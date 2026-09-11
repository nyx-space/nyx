use anise::constants::frames::EARTH_J2000;
use anise::prelude::Almanac;
use hifitime::prelude::*;
use indexmap::IndexMap;
use nyx_space::io::ExportCfg;
use nyx_space::md::prelude::*;
use nyx_space::od::DopplerConfig;
use nyx_space::od::msr::{IntegrationRef, Measurement, MeasurementType, TrackingDataArc};
use nyx_space::od::prelude::*;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::str::FromStr;
use std::sync::Arc;

#[test]
fn test_measurement_struct() {
    let epoch = Epoch::from_gregorian_utc_at_midnight(2025, 1, 1);
    let msr = Measurement::new("DSS-14".to_string(), epoch);
    assert_eq!(msr.tracker, "DSS-14");
    assert_eq!(msr.epoch, epoch);
    assert!(msr.doppler_config.is_none());

    // with() with Doppler and None should default DopplerConfig
    let msr_doppler = msr.clone().with(MeasurementType::Doppler, 1.234, None);
    assert_eq!(msr_doppler.doppler_config, Some(DopplerConfig::default()));

    // with() with custom DopplerConfig
    let custom_cfg = DopplerConfig {
        integration_time: 10 * Unit::Second,
        integration_ref: IntegrationRef::Start,
    };
    let msr_custom = msr
        .clone()
        .with(MeasurementType::Doppler, 1.234, Some(custom_cfg));
    assert_eq!(msr_custom.doppler_config, Some(custom_cfg));

    // with_doppler_config
    let msr_with_cfg = msr.with_doppler_config(Some(custom_cfg));
    assert_eq!(msr_with_cfg.doppler_config, Some(custom_cfg));
}

#[test]
fn test_tdm_integration_time_and_ref_reading() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  INTEGRATION_TIME = 15.0
  INTEGRATION_REF = START
META_STOP
DATA_START
  DOPPLER_INTEGRATED = 2025-01-01T00:00:00.000 0.123456
DATA_STOP
"#;

    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_tdm_integr_time_ref.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    File::create(&path)
        .unwrap()
        .write_all(tdm_content.as_bytes())
        .unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    assert_eq!(
        msr.doppler_config,
        Some(DopplerConfig {
            integration_time: 15 * Unit::Second,
            integration_ref: IntegrationRef::Start,
        })
    );
}

#[test]
fn test_tdm_integration_interval_reading() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  INTEGRATION_INTERVAL = 30.0
  INTEGRATION_REF = END
META_STOP
DATA_START
  DOPPLER_INTEGRATED = 2025-01-01T00:00:00.000 0.123456
DATA_STOP
"#;

    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_tdm_integr_interval.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    File::create(&path)
        .unwrap()
        .write_all(tdm_content.as_bytes())
        .unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    assert_eq!(
        msr.doppler_config,
        Some(DopplerConfig {
            integration_time: 30 * Unit::Second,
            integration_ref: IntegrationRef::End,
        })
    );
}

#[test]
fn test_tdm_integration_time_duration_string() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  INTEGRATION_TIME = 2 min
  INTEGRATION_REF = MIDDLE
META_STOP
DATA_START
  DOPPLER_INTEGRATED = 2025-01-01T00:00:00.000 0.123456
DATA_STOP
"#;

    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_tdm_integr_string.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    File::create(&path)
        .unwrap()
        .write_all(tdm_content.as_bytes())
        .unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    assert_eq!(
        msr.doppler_config,
        Some(DopplerConfig {
            integration_time: 2 * Unit::Minute,
            integration_ref: IntegrationRef::Middle,
        })
    );
}

#[test]
fn test_tdm_export_and_import_roundtrip() {
    let epoch = Epoch::from_gregorian_utc_at_midnight(2025, 1, 1);
    let custom_cfg = DopplerConfig {
        integration_time: 45 * Unit::Second,
        integration_ref: IntegrationRef::Start,
    };

    let mut data = IndexMap::new();
    data.insert(MeasurementType::Doppler, 0.5);

    let msr = Measurement {
        tracker: "DSS14".to_string(),
        epoch,
        data,
        rejected: false,
        doppler_config: Some(custom_cfg),
    };

    let arc = TrackingDataArc {
        measurements: vec![msr],
        source: None,
        moduli: None,
        force_reject: false,
    };

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_tdm_roundtrip.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();

    arc.to_tdm_file(
        "MySpacecraft".to_string(),
        None,
        &path,
        ExportCfg::default(),
    )
    .unwrap();

    let arc_read = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc_read.len(), 1);
    let msr_read = arc_read.measurements.first().unwrap();

    assert_eq!(msr_read.doppler_config, Some(custom_cfg));
}

#[test]
fn test_parquet_export_and_import_roundtrip() {
    let epoch = Epoch::from_gregorian_utc_at_midnight(2025, 1, 1);
    let custom_cfg = DopplerConfig {
        integration_time: 25 * Unit::Second,
        integration_ref: IntegrationRef::End,
    };

    let mut data = IndexMap::new();
    data.insert(MeasurementType::Range, 10000.0);
    data.insert(MeasurementType::Doppler, 0.5);

    let msr = Measurement {
        tracker: "DSS14".to_string(),
        epoch,
        data,
        rejected: false,
        doppler_config: Some(custom_cfg),
    };

    let arc = TrackingDataArc {
        measurements: vec![msr],
        source: None,
        moduli: None,
        force_reject: false,
    };

    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_parquet_roundtrip.parquet");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();

    arc.to_parquet_simple(&path).unwrap();

    let arc_read = TrackingDataArc::from_parquet(&path).unwrap();
    assert_eq!(arc_read.len(), 1);
    let msr_read = arc_read.measurements.first().unwrap();

    assert_eq!(msr_read.tracker, "DSS14");
    assert_eq!(msr_read.epoch, epoch);
    assert_eq!(msr_read.doppler_config, Some(custom_cfg));
    assert_eq!(msr_read.data.get(&MeasurementType::Range), Some(&10000.0));
    assert_eq!(msr_read.data.get(&MeasurementType::Doppler), Some(&0.5));
}

fn load_test_almanac() -> Arc<Almanac> {
    let manifest_dir: PathBuf = [env!("CARGO_MANIFEST_DIR"), "../data/01_planetary"]
        .iter()
        .collect();

    let almanac = Almanac::new(&manifest_dir.join("pck08.pca").to_string_lossy())
        .unwrap()
        .load(
            &manifest_dir
                .join("earth_longterm_000101_251211_250915.bpc")
                .to_string_lossy(),
        )
        .unwrap()
        .load(
            &manifest_dir
                .join("earth_latest_high_prec.bpc")
                .to_string_lossy(),
        )
        .unwrap()
        .load(&manifest_dir.join("de440s.bsp").to_string_lossy())
        .unwrap();

    Arc::new(almanac)
}

#[test]
fn test_ground_station_measure_uses_measurement_doppler_config() {
    let almanac = load_test_almanac();
    let earth_frame = almanac.frame_info(EARTH_J2000).unwrap();

    let epoch = Epoch::from_str("2023-02-22T19:18:17.16 UTC").unwrap();
    let orbit =
        Orbit::try_keplerian_altitude(500.0, 1e-3, 30.0, 45.0, 75.0, 23.4, epoch, earth_frame)
            .unwrap();

    let (_, traj) = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()))
        .with(Spacecraft::builder().orbit(orbit).build(), almanac.clone())
        .for_duration_with_traj(1.hours())
        .unwrap();

    let mut gs = GroundStation::from_point("TestGS".to_string(), 0.0, 0.0, 0.0, EARTH_J2000);
    let doppler_cfg = DopplerConfig {
        integration_time: 10 * Unit::Second,
        integration_ref: IntegrationRef::Middle,
    };
    gs = gs
        .with_doppler_config(Some(doppler_cfg))
        .with_msr_type(MeasurementType::Doppler, StochasticNoise::default());

    let msr_opt = gs
        .measure(epoch + 10.minutes(), &traj, None, &almanac)
        .unwrap();
    if let Some(msr) = msr_opt {
        assert_eq!(msr.doppler_config, Some(doppler_cfg));
        assert!(msr.data.contains_key(&MeasurementType::Doppler));
    }
}
