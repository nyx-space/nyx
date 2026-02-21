use nyx_space::od::simulator::{Cadence, Handoff, Scheduler, Strand, TrkConfig};
use hifitime::{Epoch, TimeUnits};
use der::{Encode, Decode};

#[test]
fn test_handoff_asn1() {
    let h = Handoff::Greedy;
    let mut buf = Vec::new();
    h.encode_to_vec(&mut buf).unwrap();
    let h2 = Handoff::from_der(&buf).unwrap();
    assert_eq!(h, h2);
}

#[test]
fn test_cadence_asn1() {
    let c = Cadence::Intermittent { on: 1.0.hours(), off: 0.5.hours() };
    let mut buf = Vec::new();
    c.encode_to_vec(&mut buf).unwrap();
    let c2 = Cadence::from_der(&buf).unwrap();
    assert_eq!(c, c2);

    let c = Cadence::Continuous;
    let mut buf = Vec::new();
    c.encode_to_vec(&mut buf).unwrap();
    let c2 = Cadence::from_der(&buf).unwrap();
    assert_eq!(c, c2);
}

#[test]
fn test_scheduler_asn1() {
    let s = Scheduler::builder()
        .handoff(Handoff::Overlap)
        .cadence(Cadence::Intermittent { on: 10.0.minutes(), off: 5.0.minutes() })
        .min_samples(5)
        .sample_alignment(1.0.seconds())
        .build();

    let mut buf = Vec::new();
    s.encode_to_vec(&mut buf).unwrap();
    let s2 = Scheduler::from_der(&buf).unwrap();
    assert_eq!(s, s2);
}

#[test]
fn test_strand_asn1() {
    let epoch = Epoch::from_gregorian_utc(2023, 1, 1, 0, 0, 0, 0);
    let s = Strand {
        start: epoch,
        end: epoch + 1.0.hours(),
    };

    let mut buf = Vec::new();
    s.encode_to_vec(&mut buf).unwrap();
    let s2 = Strand::from_der(&buf).unwrap();

    assert_eq!(s, s2);
}

#[test]
fn test_trkconfig_asn1() {
    let epoch = Epoch::from_gregorian_utc(2023, 1, 1, 0, 0, 0, 0);
    let strand = Strand {
        start: epoch,
        end: epoch + 1.0.hours(),
    };

    let cfg = TrkConfig::builder()
        .sampling(10.0.seconds())
        .strands(vec![strand])
        .build();

    let mut buf = Vec::new();
    cfg.encode_to_vec(&mut buf).unwrap();
    let cfg2 = TrkConfig::from_der(&buf).unwrap();

    assert_eq!(cfg, cfg2);
}
