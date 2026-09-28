use anise::constants::SPEED_OF_LIGHT_KM_S;
use hifitime::prelude::*;
use indexmap::IndexMap;
use nyx_space::io::ExportCfg;
use nyx_space::od::DopplerConfig;
use nyx_space::od::msr::{IntegrationRef, Measurement, MeasurementType, TrackingDataArc};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[test]
fn test_multi_section_multiple_trackers() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
COMMENT First pass from DSS14
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
META_STOP
DATA_START
  RANGE = 2023-01-01T00:00:00.000 10000.0
  RANGE = 2023-01-01T00:01:00.000 10020.0
DATA_STOP

COMMENT Second pass from DSS65
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS65
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
META_STOP
DATA_START
  RANGE = 2023-01-01T02:00:00.000 20000.0
  RANGE = 2023-01-01T02:01:00.000 20040.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_multi_trackers.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 4);

    let dss14_msrs = arc.clone().filter_by_tracker("DSS14".to_string());
    assert_eq!(dss14_msrs.len(), 2);
    assert_eq!(
        *dss14_msrs.measurements[0]
            .data
            .get(&MeasurementType::Range)
            .unwrap(),
        5000.0
    );
    assert_eq!(
        *dss14_msrs.measurements[1]
            .data
            .get(&MeasurementType::Range)
            .unwrap(),
        5010.0
    );

    let dss65_msrs = arc.filter_by_tracker("DSS65".to_string());
    assert_eq!(dss65_msrs.len(), 2);
    assert_eq!(
        *dss65_msrs.measurements[0]
            .data
            .get(&MeasurementType::Range)
            .unwrap(),
        10000.0
    );
    assert_eq!(
        *dss65_msrs.measurements[1]
            .data
            .get(&MeasurementType::Range)
            .unwrap(),
        10020.0
    );
}

#[test]
fn test_multi_section_two_way_and_one_way_concurrent() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
META_STOP
DATA_START
  RANGE                 = 2023-01-01T00:00:00.000 10000.0
  DOPPLER_INSTANTANEOUS = 2023-01-01T00:00:00.000 0.2
DATA_STOP

META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2
  ANGLE_TYPE = AZEL
META_STOP
DATA_START
  ANGLE_1 = 2023-01-01T00:00:00.000 45.0
  ANGLE_2 = 2023-01-01T00:00:00.000 60.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_two_way_one_way.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    // Because both sections share the same tracker and epoch, they should be merged into 1 measurement record
    assert_eq!(arc.len(), 1);

    let msr = arc.measurements.first().unwrap();
    assert_eq!(msr.tracker, "DSS14");
    assert_eq!(msr.epoch, Epoch::from_gregorian_utc_at_midnight(2023, 1, 1));
    assert_eq!(*msr.data.get(&MeasurementType::Range).unwrap(), 5000.0);
    assert_eq!(*msr.data.get(&MeasurementType::Doppler).unwrap(), 0.1);
    assert_eq!(*msr.data.get(&MeasurementType::Azimuth).unwrap(), 45.0);
    assert_eq!(*msr.data.get(&MeasurementType::Elevation).unwrap(), 60.0);
}

#[test]
fn test_multi_section_different_timescales() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
META_STOP
DATA_START
  RANGE = 2023-01-01T00:00:00.000 10000.0
DATA_STOP

META_START
  TIME_SYSTEM = TAI
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
META_STOP
DATA_START
  RANGE = 2023-01-01T00:00:00.000 12000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_timescales.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 2);

    let utc_epoch = Epoch::from_gregorian_utc_at_midnight(2023, 1, 1);
    let tai_epoch = Epoch::from_gregorian_tai_at_midnight(2023, 1, 1);

    // TAI is ahead of UTC, so tai_epoch occurs earlier in physical time
    assert_eq!(arc.measurements[0].epoch, tai_epoch);
    assert_eq!(arc.measurements[1].epoch, utc_epoch);
}

#[test]
fn test_multi_section_frequency_ramps() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  TURNAROUND_NUMERATOR = 240
  TURNAROUND_DENOMINATOR = 221
META_STOP
DATA_START
  TRANSMIT_FREQ      = 2023-02-22T19:18:17.160 2100000000.0
  TRANSMIT_FREQ_RATE = 2023-02-22T19:18:17.160 1.0
  RECEIVE_FREQ       = 2023-02-22T19:18:27.160 2280541478.587843
DATA_STOP

META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS65
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  TURNAROUND_NUMERATOR = 240
  TURNAROUND_DENOMINATOR = 221
META_STOP
DATA_START
  TRANSMIT_FREQ      = 2023-02-22T20:18:17.160 2100000000.0
  TRANSMIT_FREQ_RATE = 2023-02-22T20:18:17.160 1.0
  RECEIVE_FREQ       = 2023-02-22T20:18:27.160 2280541478.587843
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_multi_ramps.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 2);
    assert!(
        arc.measurements[0]
            .data
            .contains_key(&MeasurementType::Doppler)
    );
    assert!(
        arc.measurements[1]
            .data
            .contains_key(&MeasurementType::Doppler)
    );
    assert_eq!(arc.measurements[0].tracker, "DSS14");
    assert_eq!(arc.measurements[1].tracker, "DSS65");
}

#[test]
fn test_multi_section_different_range_units_and_corrections() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
  CORRECTION_RANGE = 10.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-01-01T00:00:00.000 1000.0
DATA_STOP

META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS65
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = ns
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-01-01T01:00:00.000 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_multi_corrections.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 2);

    // Section 1: (1000.0 + 10.0) / 2 = 505.0 km
    let msr1 = &arc.measurements[0];
    assert_eq!(*msr1.data.get(&MeasurementType::Range).unwrap(), 505.0);

    // Section 2: (10000.0 + 5000.0) * 1e-9 * c / 2
    let msr2 = &arc.measurements[1];
    let expected_ns = (15000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert!((*msr2.data.get(&MeasurementType::Range).unwrap() - expected_ns).abs() < 1e-10);
}

#[test]
fn test_multi_section_export_and_import_roundtrip() {
    let epoch = Epoch::from_gregorian_utc_at_midnight(2025, 1, 1);
    let doppler_cfg = DopplerConfig {
        integration_time: 10 * Unit::Second,
        integration_ref: IntegrationRef::Middle,
    };

    let mut data = IndexMap::new();
    data.insert(MeasurementType::Range, 10000.0);
    data.insert(MeasurementType::Doppler, 0.5);
    data.insert(MeasurementType::Azimuth, 30.0);
    data.insert(MeasurementType::Elevation, 45.0);

    let msr = Measurement {
        tracker: "DSS14".to_string(),
        epoch,
        data,
        rejected: false,
        doppler_config: Some(doppler_cfg),
    };

    let arc = TrackingDataArc {
        measurements: vec![msr],
        source: None,
        moduli: None,
        force_reject: false,
    };

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../target/test_multi_roundtrip.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();

    arc.to_tdm_file(
        "MySpacecraft".to_string(),
        None,
        &path,
        ExportCfg::default(),
    )
    .unwrap();

    // Verify the written file has 2 META sections (two-way and one-way)
    let content = std::fs::read_to_string(&path).unwrap();
    let meta_starts = content.matches("META_START").count();
    let data_starts = content.matches("DATA_START").count();
    assert_eq!(meta_starts, 2);
    assert_eq!(data_starts, 2);

    // Read back
    let arc_read = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc_read.len(), 1);

    let msr_read = arc_read.measurements.first().unwrap();
    assert_eq!(
        *msr_read.data.get(&MeasurementType::Range).unwrap(),
        10000.0
    );
    assert_eq!(*msr_read.data.get(&MeasurementType::Doppler).unwrap(), 0.5);
    assert_eq!(*msr_read.data.get(&MeasurementType::Azimuth).unwrap(), 30.0);
    assert_eq!(
        *msr_read.data.get(&MeasurementType::Elevation).unwrap(),
        45.0
    );
    assert_eq!(msr_read.doppler_config, Some(doppler_cfg));
}

#[test]
fn test_multi_section_syntax_errors() {
    let test_cases = [
        // Nested META_START
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
META_START
META_STOP
DATA_START
DATA_STOP"#,
            "nested META_START",
        ),
        // META_START inside DATA
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
META_STOP
DATA_START
META_START
DATA_STOP"#,
            "META_START inside DATA",
        ),
        // META_STOP without META_START
        (
            r#"CCSDS_TDM_VERS = 2.0
META_STOP
DATA_START
DATA_STOP"#,
            "META_STOP without META_START",
        ),
        // DATA_START without metadata
        (
            r#"CCSDS_TDM_VERS = 2.0
DATA_START
DATA_STOP"#,
            "DATA_START without metadata",
        ),
        // DATA_START inside metadata
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
DATA_START
META_STOP
DATA_STOP"#,
            "DATA_START inside metadata",
        ),
        // Nested DATA_START
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
META_STOP
DATA_START
DATA_START
DATA_STOP"#,
            "nested DATA_START",
        ),
        // DATA_STOP without DATA_START
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
META_STOP
DATA_STOP"#,
            "DATA_STOP without DATA_START",
        ),
        // Unterminated META_START at EOF
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC"#,
            "unterminated META_START at EOF",
        ),
        // Unterminated DATA_START at EOF
        (
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
META_STOP
DATA_START
  RANGE = 2023-01-01T00:00:00.000 1000.0"#,
            "unterminated DATA_START at EOF",
        ),
        // Metadata keyword outside META_START
        (
            r#"CCSDS_TDM_VERS = 2.0
TIME_SYSTEM = UTC
META_START
META_STOP
DATA_START
DATA_STOP"#,
            "metadata keyword outside META_START",
        ),
    ];

    for (i, (tdm_content, label)) in test_cases.iter().enumerate() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("target/test_syntax_err_{i}.tdm"));
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let mut file = File::create(&path).unwrap();
        file.write_all(tdm_content.as_bytes()).unwrap();

        let res = TrackingDataArc::from_tdm(&path, None);
        assert!(
            res.is_err(),
            "Expected failure for test case '{label}', but got Ok"
        );
    }
}
