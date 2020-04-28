use super::rv::Distribution;
use super::serde_derive::Deserialize;
use std::collections::HashMap;

#[derive(Deserialize)]
pub struct Scenario {
    distr: HashMap<String, Distribution>,
}

#[test]
fn test_deser_scenario() {
    extern crate toml;

    let _md_scen: Scenario = toml::from_str(
        r#"[distr.std_normal]
        mean = 0.0
        std_dev = 1.0
        
        [distr.my_exp.exponential]
        lambda = 0.5
        "#,
    )
    .unwrap();
}
