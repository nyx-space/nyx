use std::fmt;

use log::info;
use nyx_space::{
    cosmic::{Bodies, Cosm},
    dynamics::{Drag, Harmonics, OrbitalDynamics, PointMasses, SolarPressure, SpacecraftDynamics},
    io::gravity::HarmonicsMem,
    md::Event,
    od::ui::GroundStation,
    propagators::Propagator,
    time::{Epoch, Unit},
    Orbit, Spacecraft,
};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scenarios {
    LunarTransfer,
    OrbitDesign,
}

impl fmt::Display for Scenarios {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), std::fmt::Error> {
        match self {
            Scenarios::LunarTransfer => write!(f, "Lunar Transfer"),
            Scenarios::OrbitDesign => write!(f, "Orbit Design"),
        }
    }
}

impl Scenarios {
    pub fn run(&self) {

        info!("running scenario: {:?}", self.to_string());
        match self {
            Scenarios::LunarTransfer => lunar_transfer(),
            Scenarios::OrbitDesign => orbit_design(),
        }
        info!("done.")
    }
}

fn lunar_transfer() {
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

    info!("loading assets");
    // Set up the harmonics first because we need to pass them to the overarching orbital dynamics
    // Load the harmonics from the JGM3 file (GMAT uses the JGM2 in this case).
    // It's gunzipped (hence `true` as the last parameter)
    let stor = match HarmonicsMem::from_cof("assets/JGM3.cof.gz", 20, 20, true) {
        Ok(stor) => stor,
        Err(error) => panic!("Problem reading harmonics file: {:?}", error),
    };
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

    info!("propagating");
    let (out, traj) = prop
        .with(sat)
        .until_event(30. * Unit::Day, &Event::periapsis())
        .unwrap();
}

fn orbit_design() {
    // Load the NASA NAIF DE438 planetary ephemeris.
    let cosm = Cosm::de438();
    // Grab the Earth Mean Equator J2000 frame
    let eme2k = cosm.frame("EME2000");
    // Set the initial start time of the scenario
    let epoch = Epoch::from_gregorian_tai_at_noon(2021, 2, 25);
    // Nearly circular orbit (ecc of 0.01), inclination of 49 degrees and TA at 30.0
    let orbit = Orbit::keplerian_altitude(500.0, 0.01, 49.0, 0.0, 0.0, 30.0, epoch, eme2k);

    // Define the landmark by specifying a name, a latitude, a longitude,
    // an altitude (in km) and a frame. Note that we're also "cloning"
    // the Cosm: don't worry, it's a shared object, so we're just cloning the
    // the reference to it in memory, and never loading it more than once.
    let landmark = GroundStation::from_point(
        "Eiffel Tower".to_string(),
        36.0544,
        112.1402,
        0.0,
        cosm.frame("IAU Earth"),
        cosm.clone(),
    );

    // Let's print this landmark to make sure we've created it correctly.
    info!("target: {:}", landmark);

    // Let's specify the force model to be two body dynamics
    // And use the default propagator setup: a variable step Runge-Kutta 8-9
    let setup = Propagator::default(OrbitalDynamics::two_body());

    // Use the setup to seed a propagator with the initial state we defined above.
    let mut prop = setup.with(orbit);
    // Now let's propagate and generate the trajectory so we can analyse it.
    let (final_state, traj) = match prop.for_duration_with_traj(1 * Unit::Day) {
        Ok(state) => state,
        Err(error) => panic!("Error during propagation: {:?}", error)
    };

    // Printing the state with `:o` will print its Keplerian elements
    info!("{:?}", final_state);
}
