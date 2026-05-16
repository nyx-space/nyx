use hifitime::prelude::*;
use nyx_space::od::msr::{MeasurementType, TrackingDataArc};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[test]
fn test_tdm_frequency_ramp() {
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
  TRANSMIT_FREQ_RATE = 2023-02-22T19:18:22.160 2.0
  RECEIVE_FREQ       = 2023-02-22T19:18:27.160 2280541478.587843
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_ramp.tdm");

    // Ensure directory exists
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();

    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();

    println!("{arc}");

    assert_eq!(arc.len(), 1);
    let (epoch, msr) = arc.measurements.iter().next().unwrap();

    assert_eq!(
        *epoch,
        Epoch::from_gregorian_utc(2023, 2, 22, 19, 18, 27, 160_000_000)
    );

    let doppler = msr
        .data
        .get(&MeasurementType::Doppler)
        .expect("Doppler measurement missing");

    // Expected calculation:
    // T0 = 19:18:17.160, F0 = 2.1e9, R0 = 1.0
    // T1 = 19:18:22.160, R1 = 2.0
    // F(T1) = F0 + R0 * (T1 - T0) = 2.1e9 + 1.0 * 5 = 2100000005.0
    // T2 = 19:18:27.160, FR2 = 2280541478.587843
    // F(T2) = F(T1) + R1 * (T2 - T1) = 2100000005.0 + 2.0 * 5 = 2100000015.0
    // k = 240/221
    // f_t_expected = F(T2) * k = 2100000015.0 * 240 / 221 = 2280543002.714932
    // doppler_shift = f_t_expected - FR2 = 2280543002.714932 - 2280541478.587843 = 1524.127089
    // rho_dot = (doppler_shift * c) / (2 * f_t_expected)
    // rho_dot = (1524.127089 * 299792.458) / (2 * 2280543002.714932) = 0.1001782921

    let expected_rho_dot = 0.1001782921;
    assert!(
        (doppler - expected_rho_dot).abs() < 1e-8,
        "Doppler value mismatch: got {doppler}, expected {expected_rho_dot}",
    );
}
