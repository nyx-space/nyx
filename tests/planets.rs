extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn naif_ids() {
    use nyx::celestia::{EARTH, JUPITER, MARS, MERCURY, NAIF, NEPTUNE, SATURN, URANUS, VENUS};

    assert_eq!(MERCURY::id(), 199);
    assert_eq!(MERCURY::barycenter(), 1);
    assert_eq!(MERCURY::satellite(1), 101); // Mercury does not have satellites, cf. documentation of `satellite`.

    assert_eq!(VENUS::id(), 299);
    assert_eq!(VENUS::barycenter(), 2);
    assert_eq!(VENUS::satellite(1), 201); // Venus does not have satellites, cf. documentation of `satellite`.

    assert_eq!(EARTH::id(), 399);
    assert_eq!(EARTH::barycenter(), 3);
    assert_eq!(EARTH::satellite(1), 301);

    assert_eq!(MARS::id(), 499);
    assert_eq!(MARS::barycenter(), 4);
    assert_eq!(MARS::satellite(1), 401);

    assert_eq!(JUPITER::id(), 599);
    assert_eq!(JUPITER::barycenter(), 5);
    assert_eq!(JUPITER::satellite(1), 501);

    assert_eq!(SATURN::id(), 699);
    assert_eq!(SATURN::barycenter(), 6);
    assert_eq!(SATURN::satellite(1), 601);

    assert_eq!(URANUS::id(), 799);
    assert_eq!(URANUS::barycenter(), 7);
    assert_eq!(URANUS::satellite(1), 701);

    assert_eq!(NEPTUNE::id(), 899);
    assert_eq!(NEPTUNE::barycenter(), 8);
    assert_eq!(NEPTUNE::satellite(1), 801);
}
