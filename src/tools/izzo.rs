extern crate ndarray;
extern crate ndarray_linalg;
extern crate num_complex;
extern crate num_traits;

use ndarray::*;
use ndarray_linalg::*;
use num_traits::Float;
use std::f64::consts::PI;

use crate::linalg::Vector3;

pub struct LambertSolver {
    mu: f64,
    r1: Array1<f64>,
    r2: Array1<f64>,
    tof: f64,
    m: i32,
    prograde: bool,
    low_path: bool,
    maxiter: i32,
    atol: f64,
    rtol: f64,
    full_output: bool,
}

impl LambertSolver {
    pub fn new(
        mu: f64,
        r1: Array1<f64>,
        r2: Array1<f64>,
        tof: f64,
        m: i32,
        prograde: bool,
        low_path: bool,
        maxiter: i32,
        atol: f64,
        rtol: f64,
        full_output: bool,
    ) -> LambertSolver {
        LambertSolver {
            mu,
            r1,
            r2,
            tof,
            m,
            prograde,
            low_path,
            maxiter,
            atol,
            rtol,
            full_output,
        }
    }

    fn compute_y(&self, x: f64, ll: f64) -> f64 {
        (1.0 - ll.powi(2) * (1.0 - x.powi(2))).sqrt()
    }

    fn compute_psi(&self, x: f64, y: f64, ll: f64) -> f64 {
        if -1.0 <= x && x < 1.0 {
            // Elliptic motion
            (x * y + ll * (1.0 - x.powi(2))).acos()
        } else if x > 1.0 {
            // Hyperbolic motion
            ((y - x * ll) * (x.powi(2) - 1.0).sqrt()).asinh()
        } else {
            // Parabolic motion
            0.0
        }
    }

    fn tof_equation(&self, x: f64, T0: f64, ll: f64, M: i32) -> f64 {
        let y = self.compute_y(x, ll);
        self.tof_equation_y(x, y, T0, ll, M)
    }

    fn tof_equation_y(&self, x: f64, y: f64, T0: f64, ll: f64, M: i32) -> f64 {
        let psi = self.compute_psi(x, y, ll);
        let T_ =
            (psi + M as f64 * std::f64::consts::PI).sqrt() / (1.0 - x.powi(2)).abs() - x + ll * y;
        T_ / (1.0 - x.powi(2)) - T0
    }

    fn tof_equation_p(&self, x: f64, y: f64, T: f64, ll: f64) -> f64 {
        (3.0 * T * x - 2.0 + 2.0 * ll.powi(3) * x / y) / (1.0 - x.powi(2))
    }

    fn tof_equation_p2(&self, x: f64, y: f64, T: f64, dT: f64, ll: f64) -> f64 {
        (3.0 * T + 5.0 * x * dT + 2.0 * (1.0 - ll.powi(2)) * ll.powi(3) / y.powi(3))
            / (1.0 - x.powi(2))
    }

    fn tof_equation_p3(&self, x: f64, y: f64, T: f64, dT: f64, ddT: f64, ll: f64) -> f64 {
        (7.0 * x * ddT + 8.0 * dT - 6.0 * (1.0 - ll.powi(2)) * ll.powi(5) * x / y.powi(5))
            / (1.0 - x.powi(2))
    }

    fn compute_T_min(&self, ll: f64, M: i32, maxiter: i32, atol: f64, rtol: f64) -> (f64, f64) {
        if ll == 1.0 {
            let x_T_min = 0.0;
            let T_min = self.tof_equation(x_T_min, 0.0, ll, M);
            (x_T_min, T_min)
        } else if M == 0 {
            let x_T_min = std::f64::INFINITY;
            let T_min = 0.0;
            (x_T_min, T_min)
        } else {
            let x_i = 0.1;
            let T_i = self.tof_equation(x_i, 0.0, ll, M);
            let x_T_min = self.halley(x_i, T_i, ll, atol, rtol, maxiter);
            let T_min = self.tof_equation(x_T_min, 0.0, ll, M);
            (x_T_min, T_min)
        }
    }

    fn initial_guess(&self, T: f64, ll: f64, M: i32, low_path: bool) -> f64 {
        if M == 0 {
            let T_0 = ll.acos() + ll * (1.0 - ll.powi(2)).sqrt() + M as f64 * std::f64::consts::PI;
            let T_1 = 2.0 * (1.0 - ll.powi(3)) / 3.0;
            if T >= T_0 {
                (T_0 / T).powf(2.0 / 3.0) - 1.0
            } else if T < T_1 {
                5.0 / 2.0 * T_1 / T * (T_1 - T) / (1.0 - ll.powi(5)) + 1.0
            } else {
                (T_0 / T).powf((T_1 / T_0).log2()) - 1.0
            }
        } else {
            let x_0l = (((M as f64 * std::f64::consts::PI + std::f64::consts::PI) / (8.0 * T))
                .powf(2.0 / 3.0)
                - 1.0)
                / (((M as f64 * std::f64::consts::PI + std::f64::consts::PI) / (8.0 * T))
                    .powf(2.0 / 3.0)
                    + 1.0);
            let x_0r = (((8.0 * T) / (M as f64 * std::f64::consts::PI)).powf(2.0 / 3.0) - 1.0)
                / (((8.0 * T) / (M as f64 * std::f64::consts::PI)).powf(2.0 / 3.0) + 1.0);
            if low_path {
                x_0l.max(x_0r)
            } else {
                x_0l.min(x_0r)
            }
        }
    }

    fn halley(&self, p0: f64, T0: f64, ll: f64, atol: f64, rtol: f64, maxiter: i32) -> f64 {
        let mut p0 = p0;
        for _ in 1..=maxiter {
            let y = self.compute_y(p0, ll);
            let fder = self.tof_equation_p(p0, y, T0, ll);
            let fder2 = self.tof_equation_p2(p0, y, T0, fder, ll);
            if fder2 == 0.0 {
                panic!("Derivative was zero");
            }
            let fder3 = self.tof_equation_p3(p0, y, T0, fder, fder2, ll);
            let p = p0 - 2.0 * fder * fder2 / (2.0 * fder2.powi(2) - fder * fder3);
            if (p - p0).abs() < rtol * p0.abs() + atol {
                return p;
            }
            p0 = p;
        }
        panic!("Failed to converge");
    }

    fn householder(
        &self,
        p0: f64,
        T0: f64,
        ll: f64,
        M: i32,
        atol: f64,
        rtol: f64,
        maxiter: i32,
    ) -> (f64, i32, f64) {
        let tic = std::time::Instant::now();
        let mut p0 = p0;
        for numiter in 1..=maxiter {
            let y = self.compute_y(p0, ll);
            let fval = self.tof_equation_y(p0, y, T0, ll, M);
            let T = fval + T0;
            let fder = self.tof_equation_p(p0, y, T, ll);
            let fder2 = self.tof_equation_p2(p0, y, T, fder, ll);
            let fder3 = self.tof_equation_p3(p0, y, T, fder, fder2, ll);
            let p = p0
                - fval
                    * ((fder.powi(2) - fval * fder2 / 2.0)
                        / (fder * (fder.powi(2) - fval * fder2) + fder3 * fval.powi(2) / 6.0));
            if (p - p0).abs() < rtol * p0.abs() + atol {
                let tac = tic.elapsed();
                let tpi = tac.as_secs_f64() / numiter as f64;
                return (p, numiter, tpi);
            }
            p0 = p;
        }
        panic!("Failed to converge");
    }

    fn izzo2015(
        &self,
        mu: f64,
        r1: Vector3<f64>,
        r2: Vector3<f64>,
        tof: f64,
        M: i32,
        prograde: bool,
        low_path: bool,
        maxiter: i32,
        atol: f64,
        rtol: f64,
        full_output: bool,
    ) -> (Vector3<f64>, Vector3<f64>, i32, f64) {
        // self.assert_parameters_are_valid(mu, r1, r2, tof, M);

        let c = r2 - r1;
        let c_norm = c.norm();
        let r1_norm = r1.norm();
        let r2_norm = r2.norm();

        let s = (r1_norm + r2_norm + c_norm) * 0.5;

        let i_r1 = r1 / r1_norm;
        let i_r2 = r2 / r2_norm;
        let i_h = i_r1.cross(i_r2);
        let i_h = i_h / i_h.norm();

        let ll = (1.0 - c_norm / s).sqrt();

        let (ll, i_t1, i_t2) = if i_h.z < 0.0 {
            ll = -ll;
            let i_t1 = i_r1.cross(i_h);
            let i_t2 = i_r2.cross(i_h);
            (ll, i_t1, i_t2)
        } else {
            let i_t1 = i_h.cross(i_r1);
            let i_t2 = i_h.cross(i_r2);
            (ll, i_t1, i_t2)
        };

        let (ll, i_t1, i_t2) = if !prograde {
            (-ll, -i_t1, -i_t2)
        } else {
            (ll, i_t1, i_t2)
        };

        let T = (2.0 * mu / s.powi(3)).sqrt() * tof;

        let (x, y, numiter, tpi) = self.find_xy(ll, T, M, maxiter, atol, rtol, low_path);

        let gamma = (mu * s / 2.0).sqrt();
        let rho = (r1_norm - r2_norm) / c_norm;
        let sigma = (1.0 - rho.powi(2)).sqrt();

        let [V_r1, V_r2, V_t1, V_t2] =
            self.reconstruct(x, y, r1_norm, r2_norm, ll, gamma, rho, sigma);

        let v1 = V_r1 * (r1 / r1_norm) + V_t1 * i_t1;
        let v2 = V_r2 * (r2 / r2_norm) + V_t2 * i_t2;

        if full_output {
            (v1, v2, numiter, tpi)
        } else {
            (v1, v2, 0, 0.0)
        }
    }

    fn reconstruct(
        &self,
        x: f64,
        y: f64,
        r1: f64,
        r2: f64,
        ll: f64,
        gamma: f64,
        rho: f64,
        sigma: f64,
    ) -> [f64; 4] {
        let V_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1;
        let V_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2;
        let V_t1 = gamma * sigma * (y + ll * x) / r1;
        let V_t2 = gamma * sigma * (y + ll * x) / r2;
        [V_r1, V_r2, V_t1, V_t2]
    }

    fn find_xy(
        &self,
        ll: f64,
        T: f64,
        M: i32,
        maxiter: i32,
        atol: f64,
        rtol: f64,
        low_path: bool,
    ) -> (f64, f64, i32, f64) {
        assert!(ll.abs() < 1.0);

        let M_max = (T / PI).floor();
        let T_00 = ll.acos() + ll * (1.0 - ll.powi(2)).sqrt();

        let M_max = if T < T_00 + M_max * PI && M_max > 0.0 {
            let (_, T_min) = self.compute_T_min(ll, M_max as i32, maxiter, atol, rtol);
            if T < T_min {
                M_max - 1.0
            } else {
                M_max
            }
        } else {
            M_max
        };

        if M as f64 > M_max {
            panic!("No feasible solution, try lower M!");
        }

        let x_0 = self.initial_guess(T, ll, M, low_path);
        let (x, numiter, tpi) = self.householder(x_0, T, ll, M, atol, rtol, maxiter);
        let y = self.compute_y(x, ll);

        (x, y, numiter, tpi)
    }

    // fn compute_y(&self, x: f64, ll: f64) -> f64 {
    //     (1.0 - ll.powi(2) * (1.0 - x.powi(2))).sqrt()
    // }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn test_izzo2015() {
        let solver = LambertSolver::new();
        let mu = 398600.4418;
        let r1 = Vector3::new(15945.34, 0.0, 0.0);
        let r2 = Vector3::new(12214.83399, 10249.46731, 0.0);
        let tof = 76.0 * 60.0;
        let M = 0;
        let prograde = true;
        let low_path = true;
        let maxiter = 35;
        let atol = 1e-10;
        let rtol = 1e-10;
        let full_output = false;

        let (v1, v2, numiter, tpi) = solver.izzo2015(
            mu,
            r1,
            r2,
            tof,
            M,
            prograde,
            low_path,
            maxiter,
            atol,
            rtol,
            full_output,
        );

        assert_eq!(v1, Vector3::new(2.05783310, 2.91599867, 0.0));
        assert_eq!(v2, Vector3::new(-3.45156905, 0.91030119, 0.0));
        assert_eq!(numiter, 7);
        assert_eq!(tpi, 0.392);
    }
}
