use super::serde_derive::Deserialize;

#[derive(Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub enum Distribution {
    Normal { mean: f64, std_dev: f64 },
    Cauchy { median: f64, scale: f64 },
    Exponential { lambda: f64 },
    Poisson { lambda: f64 },
}

#[test]
fn test_deser_distr() {
    extern crate toml;

    let _std_norm: Distribution = toml::from_str(
        r#"[normal]
        mean = 0.0
        std_dev = 1.0"#,
    )
    .unwrap();

    let _cauchy: Distribution = toml::from_str(
        r#"[cauchy]
        median = 0.0
        scale = 1.0"#,
    )
    .unwrap();

    let _my_exp: Distribution = toml::from_str(
        r#"[exponential]
        lambda = 0.5"#,
    )
    .unwrap();

    let _fish: Distribution = toml::from_str(
        r#"[poisson]
        lambda = 10.0
    "#,
    )
    .unwrap();
}

#[test]
fn test_deser_distr_multi() {
    use std::collections::HashMap;
    extern crate toml;

    #[derive(Deserialize)]
    struct MapRv {
        pub rvs: HashMap<String, Distribution>,
    }

    let _as_map: MapRv = toml::from_str(
        r#"rvs = { one = { normal = { mean = 0.0, std_dev = 0.2 } }, two = { poisson = { lambda = 10.0 } } }"#,
    )
    .unwrap();
}
