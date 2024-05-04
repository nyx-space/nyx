mod cosmic;
mod mission_design;
mod monte_carlo;
mod orbit_determination;
mod propagation;
mod propulsion;

use anise::prelude::{Almanac, MetaAlmanac};

pub fn test_almanac() -> Almanac {
    use std::path::PathBuf;

    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")
        .or_else(Ok("."))
        .unwrap();

    MetaAlmanac::new(
        (PathBuf::from(manifest_dir).join("data/ci_almanac.dhall"))
            .to_string_lossy()
            .to_owned(),
    )
    .process()
    .unwrap()
}
