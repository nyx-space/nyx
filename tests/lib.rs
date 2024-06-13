mod cosmic;
mod mission_design;
mod monte_carlo;
mod orbit_determination;
mod propagation;
mod propulsion;

use std::sync::Arc;

use anise::{
    constants::celestial_objects::{EARTH, MOON, SUN},
    prelude::Almanac,
};
use propagation::{GMAT_EARTH_GM, GMAT_MOON_GM, GMAT_SUN_GM};

pub fn test_almanac() -> Almanac {
    use std::path::PathBuf;

    let manifest_dir =
        PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap_or(".".to_string()));

    Almanac::new(
        &manifest_dir
            .clone()
            .join("data/de440s.bsp")
            .to_string_lossy(),
    )
    .unwrap()
    .load(
        &manifest_dir
            .clone()
            .join("data/pck08.pca")
            .to_string_lossy(),
    )
    .unwrap()
    .load(
        &manifest_dir
            .join("data/earth_latest_high_prec.bpc")
            .to_string_lossy(),
    )
    .unwrap()
}

pub fn test_almanac_arcd() -> Arc<Almanac> {
    Arc::new(test_almanac())
}

pub fn test_almanac_gmat_arcd() -> Arc<Almanac> {
    let mut almanac = test_almanac();
    // Update GM values
    let mut earth = almanac.planetary_data.get_by_id(EARTH).unwrap();
    earth.mu_km3_s2 = GMAT_EARTH_GM;
    almanac.planetary_data.set_by_id(EARTH, earth).unwrap();

    let mut sun = almanac.planetary_data.get_by_id(SUN).unwrap();
    sun.mu_km3_s2 = GMAT_SUN_GM;
    almanac.planetary_data.set_by_id(SUN, sun).unwrap();

    let mut moon = almanac.planetary_data.get_by_id(MOON).unwrap();
    moon.mu_km3_s2 = GMAT_MOON_GM;
    almanac.planetary_data.set_by_id(MOON, moon).unwrap();

    Arc::new(almanac)
}
