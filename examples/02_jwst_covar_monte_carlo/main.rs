#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{JUPITER_BARYCENTER, MOON, SUN},
        frames::{EARTH_J2000, MOON_J2000},
    },
};
use hifitime::{TimeUnits, Unit};
use nyx::{
    cosmic::{eclipse::EclipseLocator, Frame, Mass, MetaAlmanac, SRPData},
    dynamics::{guidance::LocalFrame, OrbitalDynamics, SolarPressure, SpacecraftDynamics},
    io::ExportCfg,
    mc::MonteCarlo,
    od::{
        prelude::{KalmanVariant, KF},
        process::SpacecraftUncertainty,
        SpacecraftODProcess,
    },
    propagators::Propagator,
    Spacecraft, State,
};

use std::{collections::BTreeMap, error::Error, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();
    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's latest MetaAlmanac.
    // For details, refer to https://github.com/nyx-space/anise/blob/master/data/latest.dhall.

    // Download the regularly update of the James Webb Space Telescope reconstucted (or definitive) ephemeris.
    // Refer to https://naif.jpl.nasa.gov/pub/naif/JWST/kernels/spk/aareadme.txt for details.
    let mut latest_jwst_ephem = MetaFile {
        uri: "https://naif.jpl.nasa.gov/pub/naif/JWST/kernels/spk/jwst_rec.bsp".to_string(),
        crc32: None,
    };
    latest_jwst_ephem.process(true)?;

    // Load this ephem in the general Almanac we're using for this analysis.
    let almanac = Arc::new(
        MetaAlmanac::latest()
            .map_err(Box::new)?
            .load_from_metafile(latest_jwst_ephem, true)?,
    );

    // By loading this ephemeris file in the ANISE GUI or ANISE CLI, we can find the NAIF ID of the JWST
    // in the BSP. We need this ID in order to query the ephemeris.
    const JWST_NAIF_ID: i32 = -170;
    // Let's build a frame in the J2000 orientation centered on the JWST.
    const JWST_J2000: Frame = Frame::from_ephem_j2000(JWST_NAIF_ID);

    // Since the ephemeris file is updated regularly, we'll just grab the latest state in the ephem.
    let (earliest_epoch, latest_epoch) = almanac.spk_domain(JWST_NAIF_ID)?;
    println!("JWST defined from {earliest_epoch} to {latest_epoch}");
    // Fetch the state, printing it in the Earth J2000 frame.
    let jwst_orbit = almanac.transform(JWST_J2000, EARTH_J2000, latest_epoch, None)?;
    println!("{jwst_orbit:x}");

    // Build the spacecraft
    // SRP area assumed to be the full sunshield and mass if 6200.0 kg, c.f. https://webb.nasa.gov/content/about/faqs/facts.html
    // SRP Coefficient of reflectivity assumed to be that of Kapton, i.e. 2 - 0.44 = 1.56, table 1 from https://amostech.com/TechnicalPapers/2018/Poster/Bengtson.pdf
    let jwst = Spacecraft::builder()
        .orbit(jwst_orbit)
        .srp(SRPData {
            area_m2: 21.197 * 14.162,
            coeff_reflectivity: 1.56,
        })
        .mass(Mass::from_dry_mass(6200.0))
        .build();

    // Build up the spacecraft uncertainty builder.
    // We can use the spacecraft uncertainty structure to build this up.
    // We start by specifying the nominal state (as defined above), then the uncertainty in position and velocity
    // in the RIC frame. We could also specify the Cr, Cd, and mass uncertainties, but these aren't accounted for until
    // Nyx can also estimate the deviation of the spacecraft parameters.
    let jwst_uncertainty = SpacecraftUncertainty::builder()
        .nominal(jwst)
        .frame(LocalFrame::RIC)
        .x_km(0.5)
        .y_km(0.3)
        .z_km(1.5)
        .vx_km_s(1e-4)
        .vy_km_s(0.6e-3)
        .vz_km_s(3e-3)
        .build();

    println!("{jwst_uncertainty}");

    // Build the Kalman filter estimate.
    // Note that we could have used the KfEstimate structure directly (as seen throughout the OD integration tests)
    // but this approach requires quite a bit more boilerplate code.
    let jwst_estimate = jwst_uncertainty.to_estimate()?;

    // Set up the spacecraft dynamics.
    // We'll use the point masses of the Earth, Sun, Jupiter (barycenter, because it's in the DE440), and the Moon.
    // We'll also enable solar radiation pressure since the James Webb has a huge and highly reflective sun shield.

    let orbital_dyn = OrbitalDynamics::point_masses(vec![MOON, SUN, JUPITER_BARYCENTER]);
    let srp_dyn = SolarPressure::new(vec![EARTH_J2000, MOON_J2000], almanac.clone())?;

    // Finalize setting up the dynamics.
    let dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn);

    // Build the propagator set up to use for the whole analysis.
    let setup = Propagator::default(dynamics);

    // All of the analysis will use this duration.
    let prediction_duration = 6.5 * Unit::Day;

    // === Covariance mapping ===
    // For the covariance mapping / prediction, we'll use the common orbit determination approach.
    // This is done by setting up a spacecraft OD process, and predicting for the analysis duration.

    let ckf = KF::new(jwst_estimate, KalmanVariant::DeviationTracking);

    // Build the propagation instance for the OD process.
    let prop = setup.with(jwst.with_stm(), almanac.clone());
    let mut odp = SpacecraftODProcess::ckf(prop, ckf, BTreeMap::new(), None, almanac.clone());

    // Define the prediction step, i.e. how often we want to know the covariance.
    let step = 1_i64.minutes();
    // Finally, predict, and export the trajectory with covariance to a parquet file.
    let od_sol = odp.predict_for(step, prediction_duration)?;
    od_sol.to_parquet("./02_jwst_covar_map.parquet", ExportCfg::default())?;

    // === Monte Carlo framework ===
    // Nyx comes with a complete multi-threaded Monte Carlo frame. It's blazing fast.

    let my_mc = MonteCarlo::new(
        jwst, // Nominal state
        jwst_estimate.to_random_variable()?,
        "02_jwst".to_string(), // Scenario name
        None, // No specific seed specified, so one will be drawn from the computer's entropy.
    );

    let num_runs = 5_000;
    let rslts = my_mc.run_until_epoch(
        setup,
        almanac.clone(),
        jwst.epoch() + prediction_duration,
        num_runs,
    );

    assert_eq!(rslts.runs.len(), num_runs);
    // Finally, export these results, computing the eclipse percentage for all of these results.

    // For all of the resulting trajectories, we'll want to compute the percentage of penumbra and umbra.
    let eclipse_loc = EclipseLocator::cislunar(almanac.clone());
    let umbra_event = eclipse_loc.to_umbra_event();
    let penumbra_event = eclipse_loc.to_penumbra_event();

    rslts.to_parquet(
        "02_jwst_monte_carlo.parquet",
        Some(vec![&umbra_event, &penumbra_event]),
        ExportCfg::default(),
        almanac,
    )?;

    Ok(())
}
