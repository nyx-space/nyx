extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn tilde_matrix() {
    use self::na::{Matrix3, Vector3};
    use nyx::utils::tilde_matrix;
    let vec = Vector3::new(1.0, 2.0, 3.0);
    let rslt = Matrix3::new(0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0);
    assert_eq!(tilde_matrix(&vec), rslt);
}

#[test]
fn diagonality() {
    use self::na::Matrix3;
    use nyx::utils::is_diagonal;

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 0.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "lower triangular"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 1.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "symmetric but not diag"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "upper triangular"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        true,
        "diagonal"
    );
}
