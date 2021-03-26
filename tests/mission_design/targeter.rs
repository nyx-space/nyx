extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::ui::*;

#[test]
fn tgt_basic_position() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);
    // let mut xi = xi_orig;

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Try to increase SMA
    let xf_desired = Orbit::keplerian(
        8_100.0,
        0.2,
        30.0,
        60.0,
        60.0,
        180.0,
        orig_dt + target_delta_t,
        eme2k,
    );

    // Define the objectives
    let objectives = vec![
        Objective {
            parameter: StateParameter::X,
            desired_value: xf_desired.x,
            tolerance: 0.1,
        },
        Objective {
            parameter: StateParameter::Y,
            desired_value: xf_desired.y,
            tolerance: 0.1,
        },
        Objective {
            parameter: StateParameter::Z,
            desired_value: xf_desired.z,
            tolerance: 0.1,
        },
    ];

    let tgt = Targeter {
        prop: Arc::new(&setup),
        objectives,
        corrector: Corrector::Velocity,
        iterations: 19,
    };

    let correction = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("{}", correction);
}

#[test]
fn tgt_basic_sma() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);
    // let mut xi = xi_orig;

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Try to increase SMA
    let xf_desired = Orbit::keplerian(
        8_100.0,
        0.2,
        30.0,
        60.0,
        60.0,
        180.0,
        orig_dt + target_delta_t,
        eme2k,
    );

    // Define the objectives
    let objectives = vec![Objective {
        parameter: StateParameter::SMA,
        desired_value: xf_desired.sma(),
        tolerance: 0.1,
    }];

    let tgt = Targeter {
        prop: Arc::new(&setup),
        objectives,
        corrector: Corrector::Velocity,
        iterations: 5,
    };

    let correction = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("{}", correction);
}
