use nyx_space::{
    cosmic::{Bodies, Cosm},
    dynamics::{Drag, Harmonics, OrbitalDynamics, PointMasses, SolarPressure, SpacecraftDynamics},
    io::gravity::HarmonicsMem,
    md::Event,
    propagators::Propagator,
    time::{Epoch, Unit},
    Orbit, Spacecraft,
};

pub enum Scenarios {
    LunarTransfer,
}

pub fn lunar_transfer() {
    // Initialize the cosm which stores the ephemeris
    let cosm = Cosm::de438();
    // Grab the frames we'll use
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU Earth");
    // Define the epoch
    let epoch = Epoch::from_gregorian_utc(2014, 7, 22, 11, 29, 10, 811_000);
    // Define the initial orbit
    let orbit = Orbit::cartesian(
        -137380.1984338506,
        75679.87867537055,
        21487.63875187856,
        -0.2324532014235503,
        -0.4462753967758019,
        0.08561205662877103,
        epoch,
        eme2k,
    );

    // Define the spacecraft
    let sat = Spacecraft::new(orbit, 1000.0, 0.0, 1.0, 15.0, 1.7, 2.2);

    // Set up the harmonics first because we need to pass them to the overarching orbital dynamics
    // Load the harmonics from the JGM3 file (GMAT uses the JGM2 in this case).
    // It's gunzipped (hence `true` as the last parameter)
    let stor = HarmonicsMem::from_cof("assets/JGM3.cof.gz", 20, 20, true).unwrap();
    // Set up the orbital dynamics: we need to specify the models one by one here
    // because the usual functions wrap the dynamics so that they can be used in a Monte Carlo
    // setup.
    let orbital_dyn = OrbitalDynamics::new(vec![
        // Note that we are only accounting for Sun, Moon and Jupiter, in addition to the integration frame's GM
        PointMasses::new(
            &[Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
            cosm.clone(),
        ),
        // Specify that these harmonics are valid only in the IAU Earth frame. We're using the
        Harmonics::from_stor(iau_earth, stor, cosm.clone()),
    ]);

    // Set up SRP and Drag second, because we need to pass them to the overarching spacecraft dynamics
    let srp = SolarPressure::default(eme2k, cosm.clone());
    let drag = Drag::std_atm1976(cosm.clone());
    // Set up the spacecraft dynamics
    let sc_dyn = SpacecraftDynamics::from_models(orbital_dyn, vec![srp, drag]);

    // Propagate until periapse
    let prop = Propagator::default(sc_dyn);

    let (out, traj) = prop
        .with(sat)
        .until_event(0.5 * Unit::Day, &Event::periapsis())
        .unwrap();
}
