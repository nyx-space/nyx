mod cosmic;
mod mission_design;
mod monte_carlo;
mod orbit_determination;
mod propagation;
mod propulsion;

use std::sync::Arc;

use anise::prelude::Almanac;

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
