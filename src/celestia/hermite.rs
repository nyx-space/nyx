use celestia::exb::Coefficients;

/// Compute the interpolation at the requested time
pub fn interpolate(time: f64, coeffs: &Coefficients) -> (f64, f64, f64) {
    let x = hermval(time, &coeffs.x);
    let y = hermval(time, &coeffs.y);
    let z = hermval(time, &coeffs.z);
    (x, y, z)
}

/// hermval is a clone of numpy's Hermite values.
///
/// Given an x value, and a list of coefficients, this will return the function's interpolation.
pub fn hermval(x: f64, c: &Vec<f64>) -> f64 {
    if c.len() == 1 {
        c[0]
    } else if c.len() == 2 {
        c[0] + c[1] * x * 2.0
    } else {
        let mut nd = c.len();
        let mut c0 = c[nd - 2];
        let mut c1 = c[nd - 1];
        for i in 3..=c.len(){
            let tmp = c0;
            nd -= 1;
            c0 = c[c.len() - i] - c1 * (2.0 * f64::from((nd - 1) as i32));
            c1 = tmp + c1 * x * 2.0;
        }
        c0 + c1 * x * 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hermval() {
        use std::f64::EPSILON;
        let coeffs = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let xes = vec![1.0, 2.0, 3.0, 4.0];
        let sols = vec![-105.0, 591.0, 5215.0, 18759.0];
        for (i, x) in xes.iter().enumerate() {
            assert!((hermval(*x, &coeffs) - sols[i]).abs() < EPSILON);
        }
    }
}
