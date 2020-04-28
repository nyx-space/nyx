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
        r#"[distribution.normal]
        mean = 0.0
        std_dev = 1.0"#,
    )
    .unwrap();

    let _cauchy: Distribution = toml::from_str(
        r#"[distribution.cauchy]
        median = 0.0
        scale = 1.0"#,
    )
    .unwrap();

    let _my_exp: Distribution = toml::from_str(
        r#"[distribution.exponential]
        lambda = 0.5"#,
    )
    .unwrap();

    let _fish: Distribution = toml::from_str(
        r#"[distribution.poisson]
        lambda = 10.0
    "#,
    )
    .unwrap();
}
