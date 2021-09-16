extern crate nyx_space as nyx;
extern crate yaserde;

use nyx::io::ccsds::tdm::Tdm;
use yaserde::de::from_str;

#[test]
fn load_orekit_examples() {
  let xml = r#"
    <?xml version="1.0" encoding="UTF-8"?>
<tdm xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     xsi:noNamespaceSchemaLocation="fdsxml-1.0-master.xsd"
     id="CCSDS_TDM_VERS" version="1.0">
  <header>
    <COMMENT>TDM example created by yyyyy-nnnA Nav Team (NASA/JPL)</COMMENT>
    <COMMENT>StarTrek 1-way data, Ka band down</COMMENT>
    <CREATION_DATE>2005-06-09T20:15:00.000</CREATION_DATE>
    <ORIGINATOR>NASA/JPL</ORIGINATOR>
  </header>
  <body>
    <segment>
      <metadata>
        <COMMENT>This is a meta-data comment</COMMENT>
        <TRACK_ID>BLAH</TRACK_ID>
        <TIME_SYSTEM>Tai</TIME_SYSTEM>
      </metadata>
    </segment>
  </body>
</tdm>
    "#;

  // Then need to movee things around for supper on ground station work
  let tdm: Tdm = from_str(xml).unwrap();

  print!("{:?}", tdm);
}
