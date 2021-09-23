extern crate nyx_space as nyx;
extern crate yaserde;

use nyx::io::ccsds::tdm::{Observation, Tdm};
use std::convert::TryFrom;
use yaserde::de::from_str;

#[test]
fn load_orekit_examples() {
    use std::fs;

    // Then need to movee things around for supper on ground station work
    let tdm: Tdm =
        from_str(&fs::read_to_string("./data/tests/ccsds/tdm/orekit_TDMExample4.xml").unwrap())
            .unwrap();

    // Check a few things from the header
    assert!(!tdm.header.comments().is_empty());
    assert_eq!(
        tdm.body.segment[0].metadata.participant(1).unwrap().name,
        "DSS-24".to_owned()
    );
    assert_eq!(
        tdm.body.segment[0].metadata.participant(2).unwrap().name,
        "yyyy-nnnA".to_owned()
    );
    assert!(
        (tdm.body.segment[0].metadata.correction_range.unwrap() - 46.7741).abs()
            < std::f64::EPSILON
    );

    // Iterate through the measurements and convert them on the fly
    for observation in &tdm.body.segment[0].data.observations {
        let obs = Observation::try_from(observation).unwrap();
        println!("{:?}", obs);
    }

    // ALl keywords demo
    let tdm: Tdm = from_str(
        &fs::read_to_string("./data/tests/ccsds/tdm/orekit_TDMExampleAllKeywordsSingleDiff.xml")
            .unwrap(),
    )
    .unwrap();

    // Iterate through the measurements and convert them on the fly
    for observation in &tdm.body.segment[0].data.observations {
        let obs = Observation::try_from(observation).unwrap();
        println!("{:?}", obs);
    }
}
