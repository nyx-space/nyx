extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::MOON_J2000;
use nyx::cosmic::{try_achieve_b_plane, BPlane, BPlaneTarget, Orbit};
use nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use nyx::md::Event;
use nyx::propagators::Propagator;
use nyx::time::Epoch;

use std::str::FromStr;
use std::sync::Arc;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

use crate::propagation::GMAT_EARTH_GM;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn val_b_plane_gmat(almanac: Arc<Almanac>) {
    // This is a reproduction of the B-plane computation from the `Ex_LunarTransfer.script` file from GMAT
    // Grab the frame
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
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
    // Propagate until periapse
    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        vec![MOON, SUN, JUPITER_BARYCENTER],
    )));

    let (out, traj) = prop
        .with(orbit.into(), almanac.clone())
        .until_event(1.1 * orbit.period().unwrap(), &Event::periapsis())
        .unwrap();

    println!("{out}\n{out:x}");

    // Here is the GMAT B Plane data from that initial state until periapse
    // We accept an error of less than 500 meters in the B Plane computation.
    struct GmatData {
        epoch: Epoch,
        b_r: f64,
        b_t: f64,
        b_angle: f64,
        b_mag: f64,
        c3: f64,
    }

    let datum = vec![
        GmatData {
            epoch: Epoch::from_str("2014-07-22 11:29:45.811000000 TAI").unwrap(),
            b_r: 103582.8265522861,
            b_t: 390899.1533788401,
            b_angle: 14.84148076633666,
            b_mag: 404390.3436889349,
            c3: -4.775126658014266,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 11:30:45.811000000 TAI").unwrap(),
            b_r: 103579.7453753034,
            b_t: 390908.5559368245,
            b_angle: 14.84071755570316,
            b_mag: 404398.6433664511,
            c3: -4.775126979358285,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 11:33:42.435000000 TAI").unwrap(),
            b_r: 103570.822948689,
            b_t: 390935.9222459147,
            b_angle: 14.8385026125719,
            b_mag: 404422.8117558807,
            c3: -4.775127922707009,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 11:43:30.459000000 TAI").unwrap(),
            b_r: 103542.6799578087,
            b_t: 391023.7382534374,
            b_angle: 14.83146370497708,
            b_mag: 404500.4949941816,
            c3: -4.775131035251997,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 12:12:10.395000000 TAI").unwrap(),
            b_r: 103473.4550130072,
            b_t: 391253.2218423206,
            b_angle: 14.81367333316203,
            b_mag: 404704.6324102616,
            c3: -4.775139888778511,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 13:38:40.347000000 TAI").unwrap(),
            b_r: 103366.2782826297,
            b_t: 391738.7611333767,
            b_angle: 14.78146419868656,
            b_mag: 405146.6949887594,
            c3: -4.775164270892401,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 16:49:00.367000000 TAI").unwrap(),
            b_r: 103513.2397578877,
            b_t: 392093.1062531007,
            b_angle: 14.78876773475208,
            b_mag: 405526.8114149536,
            c3: -4.775204578014013,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 19:57:01.440000000 TAI").unwrap(),
            b_r: 103980.0635734016,
            b_t: 391969.5957286378,
            b_angle: 14.85699344872455,
            b_mag: 405526.8395512193,
            c3: -4.775224051112321,
        },
        GmatData {
            epoch: Epoch::from_str("2014-07-22 22:48:32.066000000 TAI").unwrap(),
            b_r: 104579.9942274809,
            b_t: 391732.3347895856,
            b_angle: 14.94753435539631,
            b_mag: 405451.8433948968,
            c3: -4.775221408088609,
        },
    ];

    // Iterate through the truth data
    for data in &datum {
        let eme2k_state = traj.at(data.epoch).unwrap().orbit;
        let state = almanac.transform_to(eme2k_state, MOON_J2000, None).unwrap();
        // NOTE: The transformed state is _not_ hyperbolic with de440s! Eccentricity is 0.17.
        // Compare with Cosm to understand why this state is no longer hyperbolic, the code looks to be identical.
        println!("EME2K = {eme2k_state}\nEME2K = {eme2k_state:x}");
        println!("STATE = {state}\nSTATE = {state:x}");
        assert!(
            dbg!(eme2k_state.c3_km2_s2().unwrap() - data.c3).abs() < 1e-5,
            "invalid c3 at {}",
            data.epoch
        );

        let b_plane = BPlane::new(state).unwrap();
        assert!(
            dbg!(b_plane.b_dot_r() - data.b_r).abs() < 1.0,
            "invalid b dot R at {}",
            data.epoch
        );
        assert!(
            dbg!(b_plane.b_dot_t() - data.b_t).abs() < 1.0,
            "invalid b dot T at {}",
            data.epoch
        );
        assert!(
            dbg!(b_plane.angle() - data.b_angle).abs() < 5e-4,
            "invalid b vector angle at {}",
            data.epoch
        );
        assert!(
            dbg!(b_plane.mag() - data.b_mag).abs() < 1.0,
            "invalid b vector angle at {}",
            data.epoch
        );
    }

    // Check some stuff for the first b plane
    let eme2k_state = traj.at(datum[0].epoch).unwrap().orbit;
    let state = almanac.transform_to(eme2k_state, MOON_J2000, None).unwrap();
    let b_plane = BPlane::new(state).unwrap();
    println!("{}\n{}", b_plane, b_plane.jacobian());
    println!("bt\n{}", b_plane.b_t);
    println!("br\n{}", b_plane.b_r);
    println!("ltof\n{}", b_plane.ltof_s);
}

#[rstest]
fn b_plane_davis(almanac: Arc<Almanac>) {
    // This is a simple test from Dr. Davis' IMD class at CU Boulder.

    // Hyperbolic orbit
    let orbit = Orbit::cartesian(
        546507.344255845,
        -527978.380486028,
        531109.066836708,
        -4.9220589268733,
        5.36316523097915,
        -5.22166308425181,
        Epoch::from_gregorian_utc_at_midnight(2016, 1, 1),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    );

    let bp = BPlane::new(orbit).unwrap();
    assert!((bp.b_dot_t() - 45892.323790).abs() < 1e-5, "incorrect B_T");
    assert!((bp.b_dot_r() - 10606.210428).abs() < 1e-5, "incorrect B_R");
    println!("{} km/s\n{}", orbit.vmag_km_s(), bp);

    // Check reciprocity between the gravity assist functions.
    let phi = orbit.vinf_turn_angle_deg(300.0).unwrap();
    let rp = orbit.vinf_periapsis_km(phi).unwrap();

    assert!(
        (300.0 - rp).abs() < 1e-10,
        "turn angle to rp reciprocity failed"
    );

    // Without an LTOF target, this uses the least squares approach.
    let (delta_v, achieved_b_plane) = try_achieve_b_plane(
        orbit,
        BPlaneTarget::from_bt_br(13135.7982982557, 5022.26511510685),
    )
    .unwrap();
    println!("DV (km/s) {delta_v} leads to {achieved_b_plane}");

    assert!((delta_v[0] - -0.25386251697606466).abs() < 1e-9);
    assert!((delta_v[1] - -0.18774460089778605).abs() < 1e-9);
    assert!((delta_v[2] - 0.046145009839345504).abs() < 1e-9);

    // BUG: LTOF targeting will fail.
    // let (delta_v, achieved_b_plane) = try_achieve_b_plane(
    //     orbit,
    //     BPlaneTarget::from_targets(
    //         13135.7982982557,
    //         5022.26511510685,
    //         1 * Unit::Day + 3 * Unit::Hour,
    //     ),
    // )
    // .unwrap();
}
