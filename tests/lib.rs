mod cosmic;
mod mission_design;
mod monte_carlo;
mod orbit_determination;
mod propagation;
mod propulsion;

use std::sync::Arc;

use anise::prelude::{Almanac, MetaAlmanac};

pub fn test_almanac() -> Almanac {
    use std::path::PathBuf;

    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap_or(".".to_string());

    MetaAlmanac::new(
        (PathBuf::from(manifest_dir).join("data/ci_almanac.dhall"))
            .to_string_lossy()
            .to_string(),
    )
    .unwrap()
    .process()
    .unwrap()
}

pub fn test_almanac_arcd() -> Arc<Almanac> {
    Arc::new(test_almanac())
}
