//! Core NRLMSISE-00 computation.
//!
//! Implements the mathematical framework from:
//! - Hedin (1987): MSIS-86 thermospheric model formulation
//! - Hedin (1991): Extension into middle/lower atmosphere
//! - Picone et al. (2002): NRLMSISE-00 updates
//!
//! All coefficient indices are 0-based (matching Rust arrays).

use super::Nrlmsise00Input;
use super::coefficients::*;
use core::f64::consts::PI;

const DEG_TO_RAD: f64 = PI / 180.0;
const DAY_ANGLE_RATE: f64 = 2.0 * PI / 365.25;
const HOURS_TO_RAD: f64 = PI / 12.0; // hours to radians
#[allow(non_upper_case_globals)]
const GAS_CONSTANT_J_kmolK: f64 = 831.4; // gas constant [J/(kmol·K)] adjusted for km

#[allow(non_upper_case_globals)]
/// Reference altitudes for lower thermosphere spline [km].
const SPLINE_ALTITUDES_km: [f64; 5] = [120.0, 110.0, 100.0, 90.0, 72.5];

/// Thermal diffusion coefficients (alpha) per species.
/// He, O, N2, O2, Ar, (unused), H, N
const ALPHA: [f64; 9] = [-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0];

/// Altitude limits [km] for turbopause mixing and composition corrections.
/// Above these altitudes, species are in pure diffusive equilibrium.
/// [He, O, N2, O2, Ar, (unused), H, N]
#[allow(non_upper_case_globals)]
const MIXING_ALT_LIMITS_km: [f64; 8] = [200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0];

// Molecular masses [amu].
// const MOLECULAR_MASS: [f64; 9] = [4.0, 16.0, 28.0, 32.0, 40.0, 1.0, 1.0, 14.0, 16.0];

#[allow(non_upper_case_globals)]
/// Atomic mass unit [g].
const ATOMIC_MASS_UNIT_g: f64 = 1.66e-24;

// Surface gravity and effective radius

/// Compute surface gravity [cm/s²] and effective Earth radius [km].
fn surface_gravity_and_radius(lat_deg: f64) -> (f64, f64) {
    let c2 = (2.0 * lat_deg * DEG_TO_RAD).cos();
    let gv_cm_s2 = 980.616 * (1.0 - 0.0026373 * c2);
    let reff_km = 2.0 * gv_cm_s2 / (3.085462e-6 + 2.27e-9 * c2) * 1.0e-5;
    (gv_cm_s2, reff_km)
}

/// Geopotential height [km].
fn geopotential_height_km(z: f64, zl: f64, re: f64) -> f64 {
    (z - zl) * (re + zl) / (re + z)
}

// Legendre polynomials

/// Compute associated Legendre polynomials P_n^m(x) up to degree 8, order 8.
/// Uses unnormalized convention (no Condon-Shortley phase).
/// plg\[m\]\[n\] = P_n^m(x).
fn compute_legendre(sin_lat: f64) -> [[f64; 9]; 9] {
    let mut plg = [[0.0f64; 9]; 9];
    let x = sin_lat;
    let c = (1.0 - x * x).sqrt(); // cos(lat)

    // P_0^0 = 1
    plg[0][0] = 1.0;
    // P_1^0 = x
    plg[0][1] = x;
    // P_1^1 = c
    plg[1][1] = c;

    // Sectoral: P_m^m = (2m-1) * c * P_{m-1}^{m-1}
    for m in 2..=8usize {
        plg[m][m] = (2 * m - 1) as f64 * c * plg[m - 1][m - 1];
    }

    // Diagonal+1: P_{m+1}^m = (2m+1) * x * P_m^m
    #[allow(clippy::needless_range_loop)]
    for m in 0..=7usize {
        plg[m][m + 1] = (2 * m + 1) as f64 * x * plg[m][m];
    }

    // General recurrence: (n-m)*P_n^m = (2n-1)*x*P_{n-1}^m - (n+m-1)*P_{n-2}^m
    #[allow(clippy::needless_range_loop)]
    for m in 0..=8usize {
        for n in (m + 2)..=8usize {
            plg[m][n] = ((2 * n - 1) as f64 * x * plg[m][n - 1]
                - (n + m - 1) as f64 * plg[m][n - 2])
                / (n - m) as f64;
        }
    }

    plg
}

// Ap magnetic activity functions

/// Saturation function for Ap geomagnetic index.
fn ap_saturation(a: f64, p24: f64, p25: f64) -> f64 {
    let abs_p24 = p24.abs();
    (a - 4.0) + (p25 - 1.0) * ((a - 4.0) + ((-abs_p24 * (a - 4.0)).exp() - 1.0) / abs_p24)
}

/// Exponential sum factor for Ap weighting.
fn ap_sum_factor(ex: f64) -> f64 {
    1.0 + (1.0 - ex.powi(19)) / (1.0 - ex) * ex.sqrt()
}

/// Geomagnetic activity function using 3-hour Ap history.
fn ap_geomagnetic_index(ex: f64, ap_array: &[f64; 7], p: &[f64]) -> f64 {
    // 7 elements: current + 6 historical 3-hour Ap indices (NRLMSISE-00 spec).
    let mut g0_vals = [0.0; 7];
    for (i, &a) in ap_array.iter().enumerate() {
        g0_vals[i] = ap_saturation(a, p[24], p[25]);
    }

    let sum = g0_vals[0]
        + g0_vals[1] * ex
        + g0_vals[2] * ex * ex
        + g0_vals[3] * ex.powi(3)
        + g0_vals[4] * ex.powi(4)
        + (g0_vals[5] * ex.powi(12) + g0_vals[6] * ex.powi(25)) * (1.0 - ex.powi(8)) / (1.0 - ex);

    sum / ap_sum_factor(ex)
}

// 150-term geographic/temporal variation

/// Evaluate the 150-term geographic/temporal variation function.
///
/// Returns the fractional perturbation for a 150-coefficient array `p`.
/// `sw` are the variation switches (all 1.0 for standard model).
fn geographic_variation_ratio(
    p: &[f64],
    input: &Nrlmsise00Input,
    sw: &[f64; 24],
    plg: &[[f64; 9]; 9],
) -> f64 {
    let mut t = [0.0f64; 14];

    let doy = input.day_of_year as f64;
    let df_sfu = input.f107_daily - input.f107_avg;
    let dfa_sfu = input.f107_avg - 150.0;

    // F10.7 modulation of asymmetric annual and diurnal
    let f1 = 1.0 + (p[47] * dfa_sfu + p[19] * df_sfu + p[20] * df_sfu * df_sfu) * sw[1].abs();
    let f2 = 1.0 + (p[49] * dfa_sfu + p[19] * df_sfu + p[20] * df_sfu * df_sfu) * sw[1].abs();

    // Local solar time harmonics
    let tloc_h = input.local_solar_time_hours;
    let tloc_rad = HOURS_TO_RAD * tloc_h;
    let ctloc = tloc_rad.cos();
    let stloc = tloc_rad.sin();
    let c2tloc = (2.0 * tloc_rad).cos();
    let s2tloc = (2.0 * tloc_rad).sin();
    let c3tloc = (3.0 * tloc_rad).cos();
    let s3tloc = (3.0 * tloc_rad).sin();

    // t[0]: F10.7 effect
    t[0] = p[19] * df_sfu * (1.0 + p[59] * dfa_sfu)
        + p[20] * df_sfu * df_sfu
        + p[21] * dfa_sfu
        + p[29] * dfa_sfu * dfa_sfu;

    // t[1]: Time-independent (latitude)
    t[1] = p[1] * plg[0][2]
        + p[2] * plg[0][4]
        + p[22] * plg[0][6]
        + p[14] * plg[0][2] * dfa_sfu * sw[1].abs()
        + p[26] * plg[0][1];

    // t[2]: Symmetric annual
    t[2] = p[18] * (DAY_ANGLE_RATE * (doy - p[31])).cos();

    // t[3]: Symmetric semiannual
    t[3] = (p[15] + p[16] * plg[0][2]) * (2.0 * DAY_ANGLE_RATE * (doy - p[17])).cos();

    // t[4]: Asymmetric annual
    t[4] = f1 * (p[9] * plg[0][1] + p[10] * plg[0][3]) * (DAY_ANGLE_RATE * (doy - p[13])).cos();

    // t[5]: Asymmetric semiannual
    t[5] = p[37] * plg[0][1] * (2.0 * DAY_ANGLE_RATE * (doy - p[38])).cos();

    // Pre-compute annual cosine for diurnal/semidiurnal modulation
    let cd14 = (DAY_ANGLE_RATE * (doy - p[13])).cos();

    // t[6]: Diurnal (with annual modulation)
    let t71 = p[11] * plg[1][2] * cd14 * sw[5].abs();
    let t72 = p[12] * plg[1][2] * cd14 * sw[5].abs();
    t[6] = f2
        * ((p[3] * plg[1][1] + p[4] * plg[1][3] + p[27] * plg[1][5] + t71) * ctloc
            + (p[6] * plg[1][1] + p[7] * plg[1][3] + p[28] * plg[1][5] + t72) * stloc);

    // t[7]: Semidiurnal (with annual modulation)
    let t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * sw[5].abs();
    let t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * sw[5].abs();
    t[7] = f2
        * ((p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc
            + (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc);

    // t[8]: Magnetic activity (Ap)
    if sw[9] == -1.0 {
        // 3-hour Ap mode
        let ap = &input.ap_array;
        if p[51] != 0.0 {
            let exp1 = (-10800.0 * p[51].abs()).exp();
            let exp1 = if exp1 > 0.99999 { 0.99999 } else { exp1 };
            let sg = ap_geomagnetic_index(exp1, ap, p);
            t[8] = (p[50] * plg[0][2] + p[96] * plg[0][4]) * sg
                + (p[53] * plg[1][3] + p[98] * plg[1][5]) * sg * ctloc;
        }
    } else {
        // Daily Ap mode with saturation
        let apd = input.ap_daily - 4.0;
        let p44 = p[43].abs().max(1.0e-5);
        let p45 = p[44];
        let apdf = apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44);
        t[8] = apdf
            * (p[32]
                + p[45] * plg[0][2]
                + p[34] * plg[0][4]
                + (p[100] * plg[0][1] + p[101] * plg[0][3] + p[102] * plg[0][5])
                    * cd14
                    * sw[5].abs()
                + (p[121] * plg[1][1] + p[122] * plg[1][3] + p[123] * plg[1][5])
                    * sw[7].abs()
                    * (HOURS_TO_RAD * (tloc_h - p[124])).cos());
    }

    // t[10]: Longitudinal
    if sw[11].abs() > 0.0 && input.longitude_deg > -1000.0 {
        let lon_rad = input.longitude_deg * DEG_TO_RAD;
        t[10] = (1.0 + p[80] * dfa_sfu * sw[1].abs())
            * ((p[64] * plg[1][2]
                + p[65] * plg[1][4]
                + p[66] * plg[1][6]
                + p[103] * plg[1][1]
                + p[104] * plg[1][3]
                + p[105] * plg[1][5]
                + sw[5].abs()
                    * (p[109] * plg[1][1] + p[110] * plg[1][3] + p[111] * plg[1][5])
                    * cd14)
                * lon_rad.cos()
                + (p[90] * plg[1][2]
                    + p[91] * plg[1][4]
                    + p[92] * plg[1][6]
                    + p[106] * plg[1][1]
                    + p[107] * plg[1][3]
                    + p[108] * plg[1][5]
                    + sw[5].abs()
                        * (p[112] * plg[1][1] + p[113] * plg[1][3] + p[114] * plg[1][5])
                        * cd14)
                    * lon_rad.sin());
    }

    // t[11]: UT and mixed UT/longitude
    // Eventually consider using hifitime's UT1 tables
    if sw[12].abs() > 0.0 {
        let sr_rad_s = 7.2722e-5; // Earth rotation rate [rad/s]

        // Pure UT variation
        t[11] = (1.0 + p[95] * plg[0][1])
            * (1.0 + p[81] * dfa_sfu * sw[1].abs())
            * (1.0 + p[119] * plg[0][1] * sw[5].abs() * cd14)
            * (p[68] * plg[0][1] + p[69] * plg[0][3] + p[70] * plg[0][5])
            * (sr_rad_s * (input.ut_seconds - p[71])).cos();

        // Mixed UT/longitude coupling
        let lon_rad = input.longitude_deg.to_radians();
        t[11] += sw[11].abs()
            * (p[76] * plg[2][3] + p[77] * plg[2][5] + p[78] * plg[2][7])
            * (sr_rad_s * (input.ut_seconds - p[79]) + 2.0 * lon_rad).cos()
            * (1.0 + p[137] * dfa_sfu * sw[1].abs());
    }

    // t[12]: Mixed UT/longitude/Ap (daily Ap mode)
    if sw[13].abs() > 0.0 {
        let sr_rad_s = 7.2722e-5;
        let apdf_local = {
            let apd = input.ap_daily - 4.0;
            let p44 = p[43].abs().max(1.0e-5);
            let p45 = p[44];
            apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
        };
        // Longitude/Ap coupling
        t[12] = apdf_local * sw[11].abs() * (1.0 + p[120] * plg[0][1])
            * ((p[60] * plg[1][2] + p[61] * plg[1][4] + p[62] * plg[1][6])
                * (DEG_TO_RAD * (input.longitude_deg - p[63])).cos())
            // Seasonal longitude/Ap coupling
            + apdf_local * sw[11].abs() * sw[5].abs()
                * (p[115] * plg[1][1] + p[116] * plg[1][3] + p[117] * plg[1][5])
                * cd14
                * (DEG_TO_RAD * (input.longitude_deg - p[118])).cos()
            // Pure UT/Ap coupling
            + apdf_local * sw[12].abs()
                * (p[83] * plg[0][1] + p[84] * plg[0][3] + p[85] * plg[0][5])
                * (sr_rad_s * (input.ut_seconds - p[75])).cos();
    }

    // t[13]: Terdiurnal (with annual modulation)
    t[13] = f2
        * ((p[39] * plg[3][3] + (p[93] * plg[3][4] + p[46] * plg[3][6]) * cd14 * sw[5].abs())
            * s3tloc
            + (p[40] * plg[3][3] + (p[94] * plg[3][4] + p[48] * plg[3][6]) * cd14 * sw[5].abs())
                * c3tloc);

    // Sum: p[30] + |sw[i+1]| * t[i]
    let mut result = p[30]; // constant offset (non-zero for some coefficient arrays)
    for i in 0..14 {
        result += sw[i + 1].abs() * t[i];
    }
    result
}

// Simplified variation for lower atmosphere (100-term arrays)

fn geographic_variation_lower_ratio(
    p: &[f64],
    input: &Nrlmsise00Input,
    sw: &[f64; 24],
    plg: &[[f64; 9]; 9],
    apdf: f64,
) -> f64 {
    let mut t = [0.0f64; 14];
    let doy = input.day_of_year as f64;
    let dfa = input.f107_avg - 150.0;

    // Pre-compute seasonal cosines
    let cd32 = (DAY_ANGLE_RATE * (doy - p[31])).cos();
    let cd18 = (2.0 * DAY_ANGLE_RATE * (doy - p[17])).cos();
    let cd14 = (DAY_ANGLE_RATE * (doy - p[13])).cos();
    let cd39 = (2.0 * DAY_ANGLE_RATE * (doy - p[38])).cos();

    let tloc_h = input.local_solar_time_hours;
    let tloc_rad = HOURS_TO_RAD * tloc_h;
    let ctloc = tloc_rad.cos();
    let stloc = tloc_rad.sin();
    let c2tloc = (2.0 * tloc_rad).cos();
    let s2tloc = (2.0 * tloc_rad).sin();
    let s3tloc = (3.0 * tloc_rad).sin();
    let c3tloc = (3.0 * tloc_rad).cos();

    // t[0]: F10.7
    t[0] = p[21] * dfa;

    // t[1]: Time independent (latitude)
    t[1] = p[1] * plg[0][2]
        + p[2] * plg[0][4]
        + p[22] * plg[0][6]
        + p[26] * plg[0][1]
        + p[14] * plg[0][3]
        + p[59] * plg[0][5];

    // t[2]: Symmetrical annual
    t[2] = (p[18] + p[47] * plg[0][2] + p[29] * plg[0][4]) * cd32;

    // t[3]: Symmetrical semiannual
    t[3] = (p[15] + p[16] * plg[0][2] + p[30] * plg[0][4]) * cd18;

    // t[4]: Asymmetrical annual
    t[4] = (p[9] * plg[0][1] + p[10] * plg[0][3] + p[20] * plg[0][5]) * cd14;

    // t[5]: Asymmetrical semiannual
    t[5] = p[37] * plg[0][1] * cd39;

    // t[6]: Diurnal (with annual modulation)
    if sw[7].abs() > 0.0 {
        let t71 = p[11] * plg[1][2] * cd14 * sw[5].abs();
        let t72 = p[12] * plg[1][2] * cd14 * sw[5].abs();
        t[6] = (p[3] * plg[1][1] + p[4] * plg[1][3] + t71) * ctloc
            + (p[6] * plg[1][1] + p[7] * plg[1][3] + t72) * stloc;
    }

    // t[7]: Semidiurnal (with annual modulation)
    if sw[8].abs() > 0.0 {
        let t81 = (p[23] * plg[2][3] + p[35] * plg[2][5]) * cd14 * sw[5].abs();
        let t82 = (p[33] * plg[2][3] + p[36] * plg[2][5]) * cd14 * sw[5].abs();
        t[7] = (p[5] * plg[2][2] + p[41] * plg[2][4] + t81) * c2tloc
            + (p[8] * plg[2][2] + p[42] * plg[2][4] + t82) * s2tloc;
    }

    // t[8]: Magnetic activity (uses apdf computed by geographic_variation)
    if sw[9].abs() > 0.0 {
        t[8] = apdf * (p[32] + p[45] * plg[0][2] * sw[2].abs());
    }

    // t[10]: Longitudinal (with seasonal modulation)
    if sw[10].abs() > 0.0 && sw[11].abs() > 0.0 && input.longitude_deg > -1000.0 {
        let lon_rad = input.longitude_deg * DEG_TO_RAD;
        t[10] = (1.0
            + plg[0][1]
                * (p[80] * sw[5].abs() * (DAY_ANGLE_RATE * (doy - p[81])).cos()
                    + p[85] * sw[6].abs() * (2.0 * DAY_ANGLE_RATE * (doy - p[86])).cos())
            + p[83] * sw[3].abs() * (DAY_ANGLE_RATE * (doy - p[84])).cos()
            + p[87] * sw[4].abs() * (2.0 * DAY_ANGLE_RATE * (doy - p[88])).cos())
            * ((p[64] * plg[1][2]
                + p[65] * plg[1][4]
                + p[66] * plg[1][6]
                + p[74] * plg[1][1]
                + p[75] * plg[1][3]
                + p[76] * plg[1][5])
                * lon_rad.cos()
                + (p[90] * plg[1][2]
                    + p[91] * plg[1][4]
                    + p[92] * plg[1][6]
                    + p[77] * plg[1][1]
                    + p[78] * plg[1][3]
                    + p[79] * plg[1][5])
                    * lon_rad.sin());
    }

    // t[13]: Terdiurnal
    if sw[14].abs() > 0.0 {
        t[13] = p[39] * plg[3][3] * s3tloc + p[40] * plg[3][3] * c3tloc;
    }

    let mut result = 0.0;
    for i in 0..14 {
        result += sw[i + 1].abs() * t[i];
    }
    result
}

// Cubic Hermite spline helpers

/// Natural cubic spline setup.
///
/// `n` = number of data points (= `SPLINE_ALTITUDES.len()` = 5 in practice).
/// Work array `u` is allocated on the stack as a fixed-size buffer.
fn cubic_spline_setup(x: &[f64], y: &[f64], n: usize, yp1: f64, ypn: f64, y2: &mut [f64]) {
    // SPLINE_ALTITUDES has 5 nodes; this is the maximum n this function sees.
    const MAX_SPLINE_NODES: usize = 8; // generous upper bound
    debug_assert!(n <= MAX_SPLINE_NODES);
    let mut u = [0.0f64; MAX_SPLINE_NODES];

    if yp1 > 0.99e30 {
        y2[0] = 0.0;
        u[0] = 0.0;
    } else {
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }

    for i in 1..n - 1 {
        let sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        let p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (6.0
            * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
            / (x[i + 1] - x[i - 1])
            - sig * u[i - 1])
            / p;
    }

    if ypn > 0.99e30 {
        y2[n - 1] = 0.0;
    } else {
        let un =
            (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
        y2[n - 1] = (un - 0.5 * u[n - 2]) / (0.5 * y2[n - 2] + 1.0);
    }

    for k in (0..n - 1).rev() {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }
}

/// Cubic spline interpolation.
fn cubic_spline_interpolate(xa: &[f64], ya: &[f64], y2a: &[f64], n: usize, x: f64) -> f64 {
    let mut klo = 0;
    let mut khi = n - 1;
    while khi - klo > 1 {
        let k = (khi + klo) / 2;
        if xa[k] > x {
            khi = k;
        } else {
            klo = k;
        }
    }
    let h = xa[khi] - xa[klo];
    let a = (xa[khi] - x) / h;
    let b = (x - xa[klo]) / h;
    a * ya[klo]
        + b * ya[khi]
        + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0
}

/// Cubic spline integral from xa[0] to x.
fn cubic_spline_integrate(xa: &[f64], ya: &[f64], y2a: &[f64], n: usize, x: f64) -> f64 {
    let mut yi = 0.0;
    let mut klo = 0;
    let mut khi = 1;

    while x > xa[klo] && khi < n {
        let xx = if x < xa[khi] { x } else { xa[khi] };
        let h = xa[khi] - xa[klo];
        let a = (xa[khi] - xx) / h;
        let b = (xx - xa[klo]) / h;
        let a2 = a * a;
        let b2 = b * b;
        yi += ((1.0 - a2) * ya[klo] / 2.0
            + b2 * ya[khi] / 2.0
            + ((-(1.0 + a2 * a2) / 4.0 + a2 / 2.0) * y2a[klo]
                + (b2 * b2 / 4.0 - b2 / 2.0) * y2a[khi])
                * h
                * h
                / 6.0)
            * h;

        klo += 1;
        khi += 1;
    }
    yi
}

// Density/temperature computation

struct DensityProfile {
    temp_k: f64,
    density_per_cm3: f64,
}

/// Compute temperature and density using Bates-Walker profile above ZA,
/// with spline profile below.
///
/// ZA = zn1[0] is the joining altitude between Bates-Walker and the spline.
/// zlb is the reference altitude for the Bates-Walker profile (where T = tlb).
/// ZA may differ from zlb (e.g., ZA = 123.435 km, zlb = 120 km).
///
/// Returns (temperature, density).
#[allow(clippy::too_many_arguments)]
fn density_temperature_profile(
    alt: f64,
    dlb: f64,
    tinf_k: f64,
    tlb: f64,
    xm: f64,
    alpha: f64,
    zlb: f64,
    s: f64,
    zn1: &[f64],
    tn1: &[f64],
    tgn1: &[f64],
    gsurf: f64,
    re: f64,
) -> DensityProfile {
    let n = zn1.len();
    let za = zn1[0]; // joining altitude (top of spline)

    // Compute Bates-Walker temperature at max(alt, za), referenced to zlb
    let z_bw = if alt > za { alt } else { za };
    let zg2_km = geopotential_height_km(z_bw, zlb, re);
    let tt = tinf_k - (tinf_k - tlb) * (-s * zg2_km).exp();

    // Above ZA: return Bates-Walker temperature and density
    if alt > za {
        if xm == 0.0 {
            return DensityProfile {
                temp_k: tt,
                density_per_cm3: 0.0,
            };
        }
        let glb = gsurf / (1.0 + zlb / re).powi(2);
        let gamma = xm * glb / (s * GAS_CONSTANT_J_kmolK * tinf_k);
        let expl = (-s * gamma * zg2_km).exp();
        let density_per_cm3 = dlb * (tlb / tt).powf(1.0 + alpha + gamma) * expl;
        return DensityProfile {
            temp_k: tt,
            density_per_cm3,
        };
    }

    // Below ZA: spline interpolation of 1/T through temperature nodes.
    //
    // The top node gets the Bates-Walker temperature and gradient at za.
    // Coordinate: x = geopotential_height(z, za) / zgdif, normalized to [0, 1].
    // x = 0 at top (za), x = 1 at bottom (72.5 km).

    // Temperature and gradient at ZA from Bates-Walker profile
    let ta = tt; // Bates-Walker temperature at za
    let dta = (tinf_k - ta) * s * ((re + zlb) / (re + za)).powi(2);

    // Clamp altitude to bottom of spline region
    let z = alt.max(zn1[n - 1]);

    let z1 = za; // top of spline
    let z2 = zn1[n - 1]; // bottom of spline
    let t1 = ta; // temperature at top (Bates-Walker at za)
    let t2 = tn1[n - 1]; // temperature at bottom

    let zgdif = geopotential_height_km(z2, z1, re); // < 0 (bottom below top)

    // n = SPLINE_ALTITUDES.len() = 5 (thermosphere–mesosphere junction nodes).
    const MAX_NODES: usize = 8; // generous upper bound
    debug_assert!(n <= MAX_NODES);
    let mut xs = [0.0f64; MAX_NODES];
    let mut ys = [0.0f64; MAX_NODES]; // 1/T values
    let mut y2out = [0.0f64; MAX_NODES];

    for k in 0..n {
        xs[k] = geopotential_height_km(zn1[k], z1, re) / zgdif; // normalized [0, 1]
        ys[k] = if k == 0 { 1.0 / t1 } else { 1.0 / tn1[k] };
    }

    // Spline end-derivatives in normalized coordinate:
    // d(1/T)/dx = d(1/T)/dζ × dζ/dx = d(1/T)/dζ × zgdif
    // d(1/T)/dζ = -(1/T²) × dT/dζ
    let yd1 = -dta / (t1 * t1) * zgdif; // top derivative
    let yd2 = -tgn1[1] / (t2 * t2) * zgdif * ((re + z2) / (re + z1)).powi(2); // bottom

    cubic_spline_setup(&xs, &ys, n, yd1, yd2, &mut y2out);

    let zg = geopotential_height_km(z, z1, re);
    let x = zg / zgdif; // normalized query point
    let y = cubic_spline_interpolate(&xs, &ys, &y2out, n, x);
    let temp_k = 1.0 / y; // temperature from 1/T spline

    if xm == 0.0 {
        return DensityProfile {
            temp_k,
            density_per_cm3: 0.0,
        };
    }

    // Density at ZA from Bates-Walker (propagated from zlb to za)
    let glb_zlb = gsurf / (1.0 + zlb / re).powi(2);
    let gamma = xm * glb_zlb / (s * GAS_CONSTANT_J_kmolK * tinf_k);
    let expl_za = (-s * gamma * zg2_km).exp();
    let densa = dlb * (tlb / ta).powf(1.0 + alpha + gamma) * expl_za;

    // Density below ZA: barometric integration through spline from za to alt
    let glb_za = gsurf / (1.0 + z1 / re).powi(2);
    let gamm = xm * glb_za * zgdif / GAS_CONSTANT_J_kmolK;

    let yi = cubic_spline_integrate(&xs, &ys, &y2out, n, x);
    let mut expl = gamm * yi;
    if expl > 50.0 {
        expl = 50.0;
    }

    let density_per_cm3 = densa * (t1 / temp_k).powf(1.0 + alpha) * (-expl).exp();
    DensityProfile {
        temp_k,
        density_per_cm3,
    }
}

/// Smooth logistic transition between diffusive and mixed densities.
fn mixing_transition(dd: f64, dm: f64, zhm: f64, xmm: f64, xm: f64) -> f64 {
    let a = zhm / (xmm - xm);

    if dm <= 0.0 || dd <= 0.0 {
        return if dd > 0.0 {
            dd
        } else if dm > 0.0 {
            dm
        } else {
            1.0
        };
    }

    let ylog = a * (dm / dd).ln();
    if ylog < -10.0 {
        return dd;
    }
    if ylog > 10.0 {
        return dm;
    }
    dd * (1.0 + ylog.exp()).powf(1.0 / a)
}

/// Chemistry/composition correction factor.
fn composition_correction(alt: f64, r: f64, h1: f64, zh: f64) -> f64 {
    let e = (alt - zh) / h1;
    if e > 70.0 {
        return 1.0;
    }
    if e < -70.0 {
        return r.exp();
    }
    (r / (1.0 + e.exp())).exp()
}

/// Chemistry/composition correction factor (version 2, with two heights).
fn composition_correction_dual(alt: f64, r: f64, h1: f64, zh: f64, h2: f64) -> f64 {
    let e1 = (alt - zh) / h1;
    let e2 = (alt - zh) / h2;
    if e1 > 70.0 || e2 > 70.0 {
        return 1.0;
    }
    if e1 < -70.0 && e2 < -70.0 {
        return r.exp();
    }
    let ex1 = if e1 < -70.0 {
        (-70.0_f64).exp()
    } else {
        e1.exp()
    };
    let ex2 = if e2 < -70.0 {
        (-70.0_f64).exp()
    } else {
        e2.exp()
    };
    (r / (1.0 + 0.5 * (ex1 + ex2))).exp()
}

// Main computation (GTS7/GTD7D equivalent)

/// Full NRLMSISE-00 computation.
///
/// Returns (densities\[9\], temperature_exo, temperature_alt).
/// Densities in cm⁻³: \[He, O, N2, O2, Ar, total_mass, H, N, anomO\].
pub fn compute(input: &Nrlmsise00Input, sw: &[f64; 24]) -> ([f64; 9], f64, f64) {

    let sin_lat = input.latitude_deg.to_radians().sin();
    let plg = compute_legendre(sin_lat);

    let (gsurf_cm_s2, re_km) = surface_gravity_and_radius(input.latitude_deg);

    let alt_km = input.altitude_km;

    // Joining altitude (ZA)
    // ZA is where the Bates-Walker profile joins the spline profile.
    // The C reference uses pdl[1][15] = 123.435 km (above the physical ZLB of 120 km).
    let za = CORRECTION_PARAMS[1][15];
    let mut zn1 = SPLINE_ALTITUDES_km;
    zn1[0] = za;

    // Exospheric temperature
    // Tinf variations are not important below za (simplification from C reference).
    // The asymptotic exospheric temperature in Kelvin
    let tinf_k = if alt_km > za {
        (TEMP_BOUNDARY[0]
            * TEMP_COEFFICIENTS_K[0]
            * (1.0 + sw[16] * geographic_variation_ratio(&TEMP_COEFFICIENTS_K, input, &sw, &plg)))
        .max(0.0)
    } else {
        TEMP_BOUNDARY[0] * TEMP_COEFFICIENTS_K[0]
    };

    // Temperature at lower boundary (120 km)
    let tlb_k = TEMP_BOUNDARY[1]
        * DENSITY_COEFFICIENTS[3][0]
        * (1.0 + sw[17] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[3], input, &sw, &plg));

    // Gradient parameter (g0) and normalized s
    // g0 is in Kelvin per km and is the vertical temperature gradient at the $120\text{ km}$
    // Gradient variations not important below bottom of spline (72.5 km)
    let g0_k_km = if alt_km > SPLINE_ALTITUDES_km[4] {
        TEMP_BOUNDARY[3]
            * GRADIENT_COEFFICIENTS[0]
            * (1.0 + sw[19] * geographic_variation_ratio(&GRADIENT_COEFFICIENTS, input, &sw, &plg))
    } else {
        TEMP_BOUNDARY[3] * GRADIENT_COEFFICIENTS[0]
    };

    let s_per_km = g0_k_km / (tinf_k - tlb_k);

    // Ap saturation (apdf)
    // geographic_variation computes this as a side-effect in the C code (static variable).
    // The saturation parameters p[43], p[44] are identical across all coefficient
    // arrays, so we compute it once here. geographic_variation_lower reuses this value for its
    // magnetic activity term.
    let apdf = {
        let apd = input.ap_daily - 4.0;
        let p44 = TEMP_COEFFICIENTS_K[43].abs().max(1.0e-5);
        let p45 = TEMP_COEFFICIENTS_K[44];
        apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
    };

    // Lower thermosphere temperature profile (TN1 nodes)
    // Nodes at: za(~123), 110, 100, 90, 72.5 km
    // tn1[0] is set inside density_temperature_profile from Bates-Walker at za, not from tlb directly.
    // Lower temperature node variations are only significant below 300 km.
    let mut tn1 = [0.0f64; 5];
    let mut tgn1 = [0.0f64; 2];

    tn1[0] = tlb_k; // placeholder; density_temperature_profile overrides with Bates-Walker temp at za
    if alt_km < 300.0 {
        tn1[1] = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[0],
                        input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[2] = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[1],
                        input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[3] = TEMP_BOUNDARY[7] * TEMP_NODE_COEFFICIENTS[2][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[2],
                        input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[4] = TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]
            / (1.0
                - sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[3],
                        input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tgn1[1] = TEMP_BOUNDARY[8]
            * MID_ATMO_COEFFICIENTS[8][0]
            * (1.0
                + sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &MID_ATMO_COEFFICIENTS[8],
                        input,
                        &sw,
                        &plg,
                        apdf,
                    ))
            * tn1[4]
            * tn1[4]
            / (TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]).powi(2);
    } else {
        tn1[1] = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0];
        tn1[2] = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0];
        tn1[3] = TEMP_BOUNDARY[7] * TEMP_NODE_COEFFICIENTS[2][0];
        tn1[4] = TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0];
        tgn1[1] = TEMP_BOUNDARY[8] * MID_ATMO_COEFFICIENTS[8][0] * tn1[4] * tn1[4]
            / (TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]).powi(2);
    }
    // tgn1[0] is unused: density_temperature_profile computes gradient at za from Bates-Walker

    // Species densities
    //
    // Each species follows the NRLMSISE-00 formulation (Picone et al., 2002):
    //   1. Base density at ZLB from geographic_variation function
    //   2. Diffusive equilibrium profile via Bates-Walker (density_temperature_profile)
    //   3. Below species-specific altitude limit (MIXING_ALT_LIMITS):
    //      a. Turbopause mixing: two-stage density + mixing_transition blending
    //      b. Composition correction: composition_correction/composition_correction_dual with CORRECTION_PARAMS-scaled parameters
    //   4. Some species have unconditional corrections (photodissociation, chemistry)
    //
    let mut d = [0.0f64; 9]; // He, O, N2, O2, Ar, total_mass, H, N, anomO

    let xmm = DENSITY_BOUNDARY[2][4]; // mean molecular mass at turbopause = 28.95
    let zlb_km = TEMP_BOUNDARY[5]; // lower boundary altitude = 120 km
    let zhm28 = DENSITY_BOUNDARY[2][3] * CORRECTION_PARAMS[1][5]; // turbopause mixing transition altitude parameter = 28.0

    // Turbopause height variation factor (latitude/season dependent)
    let zhf = CORRECTION_PARAMS[1][24]
        * (1.0
            + sw[5]
                * CORRECTION_PARAMS[0][24]
                * (input.latitude_deg * DEG_TO_RAD).sin()
                * (DAY_ANGLE_RATE * (input.day_of_year as f64 - TEMP_COEFFICIENTS_K[13])).cos());

    // N2 reference density at ZLB
    let g28 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[2], input, &sw, &plg);
    let db28 = DENSITY_BOUNDARY[2][0] * g28.exp() * DENSITY_COEFFICIENTS[2][0];

    // N2 mixing reference at turbopause height (shared by all species for rl)
    let zh28 = DENSITY_BOUNDARY[2][2] * zhf;
    let b28 = density_temperature_profile(
        zh28,
        db28,
        tinf_k,
        tlb_k,
        28.0 - xmm,
        ALPHA[2] - 1.0,
        zlb_km,
        s_per_km,
        &zn1,
        &tn1,
        &tgn1,
        gsurf_cm_s2,
        re_km,
    )
    .density_per_cm3;

    // He
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[0], input, &sw, &plg);
        let db04 = DENSITY_BOUNDARY[0][0] * g1.exp() * DENSITY_COEFFICIENTS[0][0];
        #[cfg(test)]
        eprintln!("  He: g1={g1:.8} db04={db04:.6e}");
        let dd = density_temperature_profile(
            alt_km,
            db04,
            tinf_k,
            tlb_k,
            4.0,
            ALPHA[0],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[0] = dd;
        if alt_km < MIXING_ALT_LIMITS_km[0] {
            // Turbopause mixing: compute mixing density via two-stage profile
            let zh04 = DENSITY_BOUNDARY[0][2];
            let b04 = density_temperature_profile(
                zh04,
                db04,
                tinf_k,
                tlb_k,
                4.0 - xmm,
                ALPHA[0] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm04 = density_temperature_profile(
                alt_km,
                b04,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[0] = mixing_transition(d[0], dm04, zhm28, xmm, 4.0);
            // Composition correction (ground mixing ratio)
            let rl = (b28 * DENSITY_BOUNDARY[0][1] / b04).ln();
            let hc04 = DENSITY_BOUNDARY[0][5] * CORRECTION_PARAMS[1][1];
            let zc04 = DENSITY_BOUNDARY[0][4] * CORRECTION_PARAMS[1][0];
            d[0] *= composition_correction(alt_km, rl, hc04, zc04);
        }
    }

    // O
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[1], input, &sw, &plg);
        let db16 = DENSITY_BOUNDARY[1][0] * g1.exp() * DENSITY_COEFFICIENTS[1][0];
        let dd = density_temperature_profile(
            alt_km,
            db16,
            tinf_k,
            tlb_k,
            16.0,
            ALPHA[1],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[1] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[1] {
            // Turbopause mixing
            let zh16 = DENSITY_BOUNDARY[1][2];
            let b16 = density_temperature_profile(
                zh16,
                db16,
                tinf_k,
                tlb_k,
                16.0 - xmm,
                ALPHA[1] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm16 = density_temperature_profile(
                alt_km,
                b16,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[1] = mixing_transition(d[1], dm16, zhm28, xmm, 16.0);
            // Diffusive equilibrium departure (F10.7-dependent)
            let rl = DENSITY_BOUNDARY[1][1]
                * CORRECTION_PARAMS[1][16]
                * (1.0 + sw[1] * CORRECTION_PARAMS[0][23] * (input.f107_avg - 150.0));
            let hc16 = DENSITY_BOUNDARY[1][5] * CORRECTION_PARAMS[1][3];
            let zc16 = DENSITY_BOUNDARY[1][4] * CORRECTION_PARAMS[1][2];
            let hc216 = DENSITY_BOUNDARY[1][5] * CORRECTION_PARAMS[1][4];
            d[1] *= composition_correction_dual(alt_km, rl, hc16, zc16, hc216);
            // Chemistry correction
            let hcc16 = DENSITY_BOUNDARY[1][7] * CORRECTION_PARAMS[1][13];
            let zcc16 = DENSITY_BOUNDARY[1][6] * CORRECTION_PARAMS[1][12];
            let rc16 = DENSITY_BOUNDARY[1][3] * CORRECTION_PARAMS[1][14];
            d[1] *= composition_correction(alt_km, rc16, hcc16, zcc16);
        }
    }

    // N2
    {
        let dd = density_temperature_profile(
            alt_km,
            db28,
            tinf_k,
            tlb_k,
            28.0,
            ALPHA[2],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[2] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[2] {
            // Turbopause mixing only (no composition correction for N2)
            let dm28_alt = density_temperature_profile(
                alt_km,
                b28,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[2] = mixing_transition(d[2], dm28_alt, zhm28, xmm, 28.0);
        }
    }

    // Temperature at altitude [K]
    let temp_alt_k = density_temperature_profile(
        alt_km,
        1.0,
        tinf_k,
        tlb_k,
        0.0,
        0.0,
        zlb_km,
        s_per_km,
        &zn1,
        &tn1,
        &tgn1,
        gsurf_cm_s2,
        re_km,
    )
    .temp_k;

    // O2
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[4], input, &sw, &plg);
        let db32 = DENSITY_BOUNDARY[3][0] * g1.exp() * DENSITY_COEFFICIENTS[4][0];
        let dd = density_temperature_profile(
            alt_km,
            db32,
            tinf_k,
            tlb_k,
            32.0,
            ALPHA[3],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[3] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[3] {
            // Turbopause mixing (O2 is heavier than mean mass; mixing_transition still applies)
            let zh32 = DENSITY_BOUNDARY[3][2];
            let b32 = density_temperature_profile(
                zh32,
                db32,
                tinf_k,
                tlb_k,
                32.0 - xmm,
                ALPHA[3] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm32 = density_temperature_profile(
                alt_km,
                b32,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[3] = mixing_transition(d[3], dm32, zhm28, xmm, 32.0);
            // Composition correction (ground mixing ratio)
            let rl = (b28 * DENSITY_BOUNDARY[3][1] / b32).ln();
            let hc32 = DENSITY_BOUNDARY[3][5] * CORRECTION_PARAMS[1][7];
            let zc32 = DENSITY_BOUNDARY[3][4] * CORRECTION_PARAMS[1][6];
            d[3] *= composition_correction(alt_km, rl, hc32, zc32);
        }
        // Photodissociation correction (all altitudes, F10.7-dependent)
        let hcc32 = DENSITY_BOUNDARY[3][7] * CORRECTION_PARAMS[1][22];
        let hcc232 = DENSITY_BOUNDARY[3][7] * CORRECTION_PARAMS[0][22];
        let zcc32 = DENSITY_BOUNDARY[3][6] * CORRECTION_PARAMS[1][21];
        let rc32 = DENSITY_BOUNDARY[3][3]
            * CORRECTION_PARAMS[1][23]
            * (1.0 + sw[1] * CORRECTION_PARAMS[0][23] * (input.f107_avg - 150.0));
        d[3] *= composition_correction_dual(alt_km, rc32, hcc32, zcc32, hcc232);
    }

    // Ar
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[5], input, &sw, &plg);
        let db40 = DENSITY_BOUNDARY[4][0] * g1.exp() * DENSITY_COEFFICIENTS[5][0];
        let dd = density_temperature_profile(
            alt_km,
            db40,
            tinf_k,
            tlb_k,
            40.0,
            ALPHA[4],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[4] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[4] {
            // Turbopause mixing
            let zh40 = DENSITY_BOUNDARY[4][2];
            let b40 = density_temperature_profile(
                zh40,
                db40,
                tinf_k,
                tlb_k,
                40.0 - xmm,
                ALPHA[4] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm40 = density_temperature_profile(
                alt_km,
                b40,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[4] = mixing_transition(d[4], dm40, zhm28, xmm, 40.0);
            // Composition correction
            let rl = (b28 * DENSITY_BOUNDARY[4][1] / b40).ln();
            let hc40 = DENSITY_BOUNDARY[4][5] * CORRECTION_PARAMS[1][9];
            let zc40 = DENSITY_BOUNDARY[4][4] * CORRECTION_PARAMS[1][8];
            d[4] *= composition_correction(alt_km, rl, hc40, zc40);
        }
    }

    // H
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[6], input, &sw, &plg);
        let db01 = DENSITY_BOUNDARY[5][0] * g1.exp() * DENSITY_COEFFICIENTS[6][0];
        let dd = density_temperature_profile(
            alt_km,
            db01,
            tinf_k,
            tlb_k,
            1.0,
            ALPHA[6],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[6] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[6] {
            // Turbopause mixing
            let zh01 = DENSITY_BOUNDARY[5][2];
            let b01 = density_temperature_profile(
                zh01,
                db01,
                tinf_k,
                tlb_k,
                1.0 - xmm,
                ALPHA[6] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm01 = density_temperature_profile(
                alt_km,
                b01,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[6] = mixing_transition(d[6], dm01, zhm28, xmm, 1.0);
            // Composition correction (ground mixing ratio)
            let rl = (b28 * DENSITY_BOUNDARY[5][1] * CORRECTION_PARAMS[1][17].abs() / b01).ln();
            let hc01 = DENSITY_BOUNDARY[5][5] * CORRECTION_PARAMS[1][11];
            let zc01 = DENSITY_BOUNDARY[5][4] * CORRECTION_PARAMS[1][10];
            d[6] *= composition_correction(alt_km, rl, hc01, zc01);
            // Chemistry correction
            let hcc01 = DENSITY_BOUNDARY[5][7] * CORRECTION_PARAMS[1][19];
            let zcc01 = DENSITY_BOUNDARY[5][6] * CORRECTION_PARAMS[1][18];
            let rc01 = DENSITY_BOUNDARY[5][3] * CORRECTION_PARAMS[1][20];
            d[6] *= composition_correction(alt_km, rc01, hcc01, zcc01);
        }
    }

    // N
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[7], input, &sw, &plg);
        let db14 = DENSITY_BOUNDARY[6][0] * g1.exp() * DENSITY_COEFFICIENTS[7][0];
        let dd = density_temperature_profile(
            alt_km,
            db14,
            tinf_k,
            tlb_k,
            14.0,
            ALPHA[7],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        d[7] = dd;
        if alt_km <= MIXING_ALT_LIMITS_km[7] {
            // Turbopause mixing
            let zh14 = DENSITY_BOUNDARY[6][2];
            let b14 = density_temperature_profile(
                zh14,
                db14,
                tinf_k,
                tlb_k,
                14.0 - xmm,
                ALPHA[7] - 1.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let dm14 = density_temperature_profile(
                alt_km,
                b14,
                tinf_k,
                tlb_k,
                xmm,
                0.0,
                zlb_km,
                s_per_km,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            d[7] = mixing_transition(d[7], dm14, zhm28, xmm, 14.0);
            // Composition correction (ground mixing ratio)
            let rl = (b28 * DENSITY_BOUNDARY[6][1] * CORRECTION_PARAMS[0][2].abs() / b14).ln();
            let hc14 = DENSITY_BOUNDARY[6][5] * CORRECTION_PARAMS[0][1];
            let zc14 = DENSITY_BOUNDARY[6][4] * CORRECTION_PARAMS[0][0];
            d[7] *= composition_correction(alt_km, rl, hc14, zc14);
            // Chemistry correction
            let hcc14 = DENSITY_BOUNDARY[6][7] * CORRECTION_PARAMS[0][4];
            let zcc14 = DENSITY_BOUNDARY[6][6] * CORRECTION_PARAMS[0][3];
            let rc14 = DENSITY_BOUNDARY[6][3] * CORRECTION_PARAMS[0][5];
            d[7] *= composition_correction(alt_km, rc14, hcc14, zcc14);
        }
    }

    // Anomalous O
    {
        let g1 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[8], input, &sw, &plg);
        let db16h = DENSITY_BOUNDARY[7][0] * g1.exp() * DENSITY_COEFFICIENTS[8][0];
        let tho = DENSITY_BOUNDARY[7][9] * CORRECTION_PARAMS[0][6];
        let dd = density_temperature_profile(
            alt_km,
            db16h,
            tho,
            tho,
            16.0,
            ALPHA[8],
            zlb_km,
            s_per_km,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        let zmho = DENSITY_BOUNDARY[7][4];
        let zsht = DENSITY_BOUNDARY[7][5];
        let zsho = GAS_CONSTANT_J_kmolK * tho / (gsurf_cm_s2 / (1.0 + zmho / re_km).powi(2) * 16.0);
        d[8] = dd * (-zsht / zsho * ((-(alt_km - zmho) / zsht).exp() - 1.0)).exp();
        if alt_km < zmho {
            d[8] = dd;
        }
    }

    // Total mass density [g/cm³] → convert to d[5]
    // All other items of `d` are in cm^-3.
    d[5] = ATOMIC_MASS_UNIT_g
        * (4.0 * d[0]       // He
            + 16.0 * d[1]   // O
            + 28.0 * d[2]   // N2
            + 32.0 * d[3]   // O2
            + 40.0 * d[4]   // Ar
            + d[6]           // H (mass=1)
            + 14.0 * d[7]   // N
            + 16.0 * d[8]); // anomalous O

    (d, tinf_k, temp_alt_k)
}

#[cfg(test)]
mod ut_nrlmsise {
    use super::*;
    use crate::dynamics::nrlmsise00::Nrlmsise00Input;

    fn test_input() -> Nrlmsise00Input {
        Nrlmsise00Input {
            day_of_year: 80,
            ut_seconds: 43200.0,
            altitude_km: 400.0,
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            local_solar_time_hours: 12.0,
            f107_daily: 150.0,
            f107_avg: 150.0,
            ap_daily: 15.0,
            ap_array: [15.0; 7],
        }
    }

    #[test]
    fn debug_compute_400km() {
        let input = test_input();
        let (d, tinf, temp) = compute(&input, &[1.0; 24]);
        eprintln!("tinf={tinf:.1} temp={temp:.1}");
        eprintln!("rho={:.4e} kg/m3", d[5] * 1000.0);
        eprintln!(
            "He={:.3e} O={:.3e} N2={:.3e} O2={:.3e}",
            d[0], d[1], d[2], d[3]
        );
        assert!(tinf > 500.0 && tinf < 2000.0, "tinf={tinf}");
        assert!(temp > 500.0, "temp={temp}");
        assert!(d[5] > 0.0, "total density must be positive");
    }

    #[test]
    fn debug_temperature_profile() {
        let input = test_input();
        let sw = [1.0f64; 24];
        let sin_lat = 0.0;
        let plg = compute_legendre(sin_lat);
        let (gsurf_cm_s2, re_km) = surface_gravity_and_radius(0.0);
        let tinf = TEMP_BOUNDARY[0]
            * TEMP_COEFFICIENTS_K[0]
            * (1.0 + sw[16] * geographic_variation_ratio(&TEMP_COEFFICIENTS_K, &input, &sw, &plg));
        let tlb = TEMP_BOUNDARY[1]
            * DENSITY_COEFFICIENTS[3][0]
            * (1.0
                + sw[17] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[3], &input, &sw, &plg));
        let g0 = TEMP_BOUNDARY[3]
            * GRADIENT_COEFFICIENTS[0]
            * (1.0
                + sw[19] * geographic_variation_ratio(&GRADIENT_COEFFICIENTS, &input, &sw, &plg));
        let s = g0 / (tinf - tlb);
        let apdf = {
            let apd = input.ap_daily - 4.0;
            let p44 = TEMP_COEFFICIENTS_K[43].abs().max(1.0e-5);
            let p45 = TEMP_COEFFICIENTS_K[44];
            apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
        };

        eprintln!("tinf={tinf:.1} tlb={tlb:.1} g0={g0:.3} s={s:.6}");

        let mut tn1 = [0.0f64; 5];
        tn1[0] = tlb;
        tn1[1] = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[0],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[2] = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[1],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[3] = TEMP_BOUNDARY[7] * TEMP_NODE_COEFFICIENTS[2][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[2],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[4] = TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]
            / (1.0
                - sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[3],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        eprintln!("tn1: {:?}", tn1);

        let mut tgn1 = [0.0f64; 2];
        tgn1[0] = g0;
        tgn1[1] = TEMP_BOUNDARY[8]
            * MID_ATMO_COEFFICIENTS[8][0]
            * (1.0
                + sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &MID_ATMO_COEFFICIENTS[8],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ))
            * tn1[4]
            * tn1[4]
            / (TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]).powi(2);
        eprintln!("tgn1: {:?}", tgn1);

        // Test temperatures at each altitude
        let zn1 = SPLINE_ALTITUDES_km;
        for alt in [
            72.5, 80.0, 90.0, 100.0, 110.0, 115.0, 120.0, 150.0, 200.0, 400.0,
        ] {
            let t = density_temperature_profile(
                alt,
                1.0,
                tinf,
                tlb,
                0.0,
                0.0,
                TEMP_BOUNDARY[5],
                s,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .temp_k;
            eprintln!("  T({alt:.1}km) = {t:.1} K");
        }
    }

    #[test]
    fn debug_species_diagnostics() {
        let sw = [1.0f64; 24];
        let sin_lat = 0.0;
        let plg = compute_legendre(sin_lat);
        let (gsurf_cm_s2, re_km) = surface_gravity_and_radius(0.0);
        let zn1 = SPLINE_ALTITUDES_km;

        let input = test_input();
        let tinf = TEMP_BOUNDARY[0]
            * TEMP_COEFFICIENTS_K[0]
            * (1.0 + sw[16] * geographic_variation_ratio(&TEMP_COEFFICIENTS_K, &input, &sw, &plg));
        let tlb = TEMP_BOUNDARY[1]
            * DENSITY_COEFFICIENTS[3][0]
            * (1.0
                + sw[17] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[3], &input, &sw, &plg));
        let g0_val = TEMP_BOUNDARY[3]
            * GRADIENT_COEFFICIENTS[0]
            * (1.0
                + sw[19] * geographic_variation_ratio(&GRADIENT_COEFFICIENTS, &input, &sw, &plg));
        let s = g0_val / (tinf - tlb);
        let zlb = TEMP_BOUNDARY[5];
        let xmm = DENSITY_BOUNDARY[2][4];
        eprintln!("tinf={tinf:.2} tlb={tlb:.2} s={s:.6} zlb={zlb} xmm={xmm}");

        let apdf = {
            let apd = input.ap_daily - 4.0;
            let p44 = TEMP_COEFFICIENTS_K[43].abs().max(1.0e-5);
            let p45 = TEMP_COEFFICIENTS_K[44];
            apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
        };

        // Temperature nodes
        let mut tn1 = [0.0f64; 5];
        let mut tgn1 = [0.0f64; 2];
        tn1[0] = tlb;
        tn1[1] = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[0],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[2] = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[1],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[3] = TEMP_BOUNDARY[7] * TEMP_NODE_COEFFICIENTS[2][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[2],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[4] = TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]
            / (1.0
                - sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[3],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tgn1[0] = g0_val;
        tgn1[1] = TEMP_BOUNDARY[8]
            * MID_ATMO_COEFFICIENTS[8][0]
            * (1.0
                + sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &MID_ATMO_COEFFICIENTS[8],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ))
            * tn1[4]
            * tn1[4]
            / (TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]).powi(2);

        // N2 base density
        let g28 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[2], &input, &sw, &plg);
        let db28 = DENSITY_BOUNDARY[2][0] * g28.exp() * DENSITY_COEFFICIENTS[2][0];
        eprintln!(
            "N2: g28={g28:.6} db28={db28:.4e} DENSITY_BOUNDARY[2][0]={} DENSITY_COEFFICIENTS[2][0]={}",
            DENSITY_BOUNDARY[2][0], DENSITY_COEFFICIENTS[2][0]
        );

        // O2 base density
        let g_o2 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[4], &input, &sw, &plg);
        let db32 = DENSITY_BOUNDARY[3][0] * g_o2.exp() * DENSITY_COEFFICIENTS[4][0];
        eprintln!(
            "O2: g_o2={g_o2:.6} db32={db32:.4e} DENSITY_BOUNDARY[3][0]={} DENSITY_COEFFICIENTS[4][0]={}",
            DENSITY_BOUNDARY[3][0], DENSITY_COEFFICIENTS[4][0]
        );

        // He base density
        let g_he = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[0], &input, &sw, &plg);
        let db04 = DENSITY_BOUNDARY[0][0] * g_he.exp() * DENSITY_COEFFICIENTS[0][0];
        eprintln!("He: g_he={g_he:.6} db04={db04:.4e}",);

        // O base density
        let g_o = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[1], &input, &sw, &plg);
        let db16 = DENSITY_BOUNDARY[1][0] * g_o.exp() * DENSITY_COEFFICIENTS[1][0];
        eprintln!("O:  g_o ={g_o:.6} db16={db16:.4e}");

        // N base density
        let g_n = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[7], &input, &sw, &plg);
        let db14 = DENSITY_BOUNDARY[6][0] * g_n.exp() * DENSITY_COEFFICIENTS[7][0];
        eprintln!("N:  g_n ={g_n:.6} db14={db14:.4e}");

        // Densities at multiple altitudes
        eprintln!("\nDensities at ZLB (120km) and above:");
        for alt in [120.0, 150.0, 200.0, 300.0, 400.0] {
            let n2 = density_temperature_profile(
                alt,
                db28,
                tinf,
                tlb,
                28.0,
                0.0,
                zlb,
                s,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let o2 = density_temperature_profile(
                alt,
                db32,
                tinf,
                tlb,
                32.0,
                0.0,
                zlb,
                s,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let he = density_temperature_profile(
                alt,
                db04,
                tinf,
                tlb,
                4.0,
                ALPHA[0],
                zlb,
                s,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .density_per_cm3;
            let t = density_temperature_profile(
                alt,
                1.0,
                tinf,
                tlb,
                0.0,
                0.0,
                zlb,
                s,
                &zn1,
                &tn1,
                &tgn1,
                gsurf_cm_s2,
                re_km,
            )
            .temp_k;
            eprintln!(
                "  {alt:.0}km: T={t:.1}K N2={n2:.4e} O2={o2:.4e} He={he:.4e} O2/N2={:.4}",
                o2 / n2
            );
        }

        // Full compute for comparison
        let (d, _, _) = compute(&input, &[1.0; 24]);
        eprintln!("\nFull compute at 400km:");
        eprintln!(
            "  He={:.4e} O={:.4e} N2={:.4e} O2={:.4e} Ar={:.4e} H={:.4e} N={:.4e}",
            d[0], d[1], d[2], d[3], d[4], d[6], d[7]
        );
        eprintln!(
            "  Total mass = {:.4e} g/cm3 = {:.4e} kg/m3",
            d[5],
            d[5] * 1000.0
        );
    }

    #[test]
    fn debug_n2_at_100km() {
        let sw = [1.0f64; 24];
        let sin_lat = 0.0;
        let plg = compute_legendre(sin_lat);
        let (gsurf_cm_s2, re_km) = surface_gravity_and_radius(0.0);
        let zn1 = SPLINE_ALTITUDES_km;

        let mut input = test_input();
        input.altitude_km = 100.0;

        let tinf = TEMP_BOUNDARY[0]
            * TEMP_COEFFICIENTS_K[0]
            * (1.0 + sw[16] * geographic_variation_ratio(&TEMP_COEFFICIENTS_K, &input, &sw, &plg));
        let tlb = TEMP_BOUNDARY[1]
            * DENSITY_COEFFICIENTS[3][0]
            * (1.0
                + sw[17] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[3], &input, &sw, &plg));
        let g0_val = TEMP_BOUNDARY[3]
            * GRADIENT_COEFFICIENTS[0]
            * (1.0
                + sw[19] * geographic_variation_ratio(&GRADIENT_COEFFICIENTS, &input, &sw, &plg));
        let s = g0_val / (tinf - tlb);
        let zlb = TEMP_BOUNDARY[5];
        let xmm = DENSITY_BOUNDARY[2][4];
        let zhm28 = DENSITY_BOUNDARY[2][4];

        let apdf = {
            let apd = input.ap_daily - 4.0;
            let p44 = TEMP_COEFFICIENTS_K[43].abs().max(1.0e-5);
            let p45 = TEMP_COEFFICIENTS_K[44];
            apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
        };

        let mut tn1 = [0.0f64; 5];
        let mut tgn1 = [0.0f64; 2];
        tn1[0] = tlb;
        tn1[1] = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[0],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[2] = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[1],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[3] = TEMP_BOUNDARY[7] * TEMP_NODE_COEFFICIENTS[2][0]
            / (1.0
                - sw[18]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[2],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tn1[4] = TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]
            / (1.0
                - sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &TEMP_NODE_COEFFICIENTS[3],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ));
        tgn1[0] = g0_val;
        tgn1[1] = TEMP_BOUNDARY[8]
            * MID_ATMO_COEFFICIENTS[8][0]
            * (1.0
                + sw[18]
                    * sw[20]
                    * geographic_variation_lower_ratio(
                        &MID_ATMO_COEFFICIENTS[8],
                        &input,
                        &sw,
                        &plg,
                        apdf,
                    ))
            * tn1[4]
            * tn1[4]
            / (TEMP_BOUNDARY[4] * TEMP_NODE_COEFFICIENTS[3][0]).powi(2);

        let alt = 100.0;
        let g28 = sw[20] * geographic_variation_ratio(&DENSITY_COEFFICIENTS[2], &input, &sw, &plg);
        let db28 = DENSITY_BOUNDARY[2][0] * g28.exp() * DENSITY_COEFFICIENTS[2][0];
        let zh28 = DENSITY_BOUNDARY[2][2];
        let b28 = density_temperature_profile(
            zh28,
            db28,
            tinf,
            tlb,
            xmm - 28.0,
            ALPHA[2] - 1.0,
            zlb,
            s,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;

        // Raw diffusive N2 at 100km
        let dd = density_temperature_profile(
            alt,
            db28,
            tinf,
            tlb,
            28.0,
            ALPHA[2],
            zlb,
            s,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;
        // Mixing density at 100km
        let dm28_alt = density_temperature_profile(
            alt,
            b28,
            tinf,
            tlb,
            xmm,
            0.0,
            zlb,
            s,
            &zn1,
            &tn1,
            &tgn1,
            gsurf_cm_s2,
            re_km,
        )
        .density_per_cm3;

        let a = zhm28 / (xmm - 28.0);
        let ratio = dm28_alt / dd;
        let result = mixing_transition(dd, dm28_alt, zhm28, xmm, 28.0);

        eprintln!("N2 at {alt}km:");
        eprintln!("  db28={db28:.4e} (at ZLB=120km)");
        eprintln!("  zh28={zh28} b28={b28:.4e}");
        eprintln!("  dd (diffusive) = {dd:.4e}");
        eprintln!("  dm28 (mixing)  = {dm28_alt:.4e}");
        eprintln!("  ratio dm/dd    = {ratio:.6}");
        eprintln!("  mixing_transition a         = {a:.4}");
        eprintln!("  mixing_transition result    = {result:.4e}");
        eprintln!("  expected (pymsis) ≈ 1.12e13 cm⁻³");
    }

    #[test]
    fn composition_correction_limits() {
        // High altitude (e >> 0): logistic → 0, exp(0) = 1.0
        let v = composition_correction(200.0, 0.5, 10.0, 50.0);
        assert!(
            (v - 1.0).abs() < 1e-6,
            "composition_correction high alt: {v}"
        );

        // Very low altitude (e < -70): should return exp(r)
        let v = composition_correction(-800.0, 0.5, 10.0, 50.0);
        assert!(
            (v - 0.5_f64.exp()).abs() < 1e-10,
            "composition_correction low alt: {v}"
        );

        // Normal range: between 1.0 and exp(r)
        let v = composition_correction(50.0, 0.5, 10.0, 50.0);
        assert!(
            v > 1.0 && v < 0.5_f64.exp(),
            "composition_correction normal: {v}"
        );
    }

    #[test]
    fn debug_f107_sensitivity() {
        let sw = [1.0f64; 24];
        let plg = compute_legendre(0.0);

        for (label, f107) in [
            ("solar_min", 70.0),
            ("solar_mod", 150.0),
            ("solar_max", 250.0),
        ] {
            let ap = if f107 < 100.0 {
                4.0
            } else if f107 < 200.0 {
                15.0
            } else {
                30.0
            };
            let input = Nrlmsise00Input {
                day_of_year: 80,
                ut_seconds: 43200.0,
                altitude_km: 100.0,
                latitude_deg: 0.0,
                longitude_deg: 0.0,
                local_solar_time_hours: 12.0,
                f107_daily: f107,
                f107_avg: f107,
                ap_daily: ap,
                ap_array: [ap; 7],
            };

            let apdf = {
                let apd = input.ap_daily - 4.0;
                let p44 = TEMP_COEFFICIENTS_K[43].abs().max(1.0e-5);
                let p45 = TEMP_COEFFICIENTS_K[44];
                apd + (p45 - 1.0) * (apd + ((-p44 * apd).exp() - 1.0) / p44)
            };

            let g7s_ptl1 = geographic_variation_lower_ratio(
                &TEMP_NODE_COEFFICIENTS[1],
                &input,
                &sw,
                &plg,
                apdf,
            );
            let tn1_2 = TEMP_BOUNDARY[2] * TEMP_NODE_COEFFICIENTS[1][0] / (1.0 - sw[18] * g7s_ptl1);

            let g7s_ptl0 = geographic_variation_lower_ratio(
                &TEMP_NODE_COEFFICIENTS[0],
                &input,
                &sw,
                &plg,
                apdf,
            );
            let tn1_1 = TEMP_BOUNDARY[6] * TEMP_NODE_COEFFICIENTS[0][0] / (1.0 - sw[18] * g7s_ptl0);

            let g7_pt = geographic_variation_ratio(&TEMP_COEFFICIENTS_K, &input, &sw, &plg);
            let tinf = TEMP_BOUNDARY[0] * TEMP_COEFFICIENTS_K[0] * (1.0 + sw[16] * g7_pt);
            let g7_pd3 = geographic_variation_ratio(&DENSITY_COEFFICIENTS[3], &input, &sw, &plg);
            let tlb = TEMP_BOUNDARY[1] * DENSITY_COEFFICIENTS[3][0] * (1.0 + sw[17] * g7_pd3);

            eprintln!("{label} (F10.7={f107}, Ap={ap}):");
            eprintln!("  apdf = {apdf:.6}");
            eprintln!(
                "  geographic_variation_lower(TEMP_NODE_COEFFICIENTS[1]) = {g7s_ptl1:.6}  → tn1[2] = {tn1_2:.2}K"
            );
            eprintln!(
                "  geographic_variation_lower(TEMP_NODE_COEFFICIENTS[0]) = {g7s_ptl0:.6}  → tn1[1] = {tn1_1:.2}K"
            );
            eprintln!("  tinf = {tinf:.1}K  tlb = {tlb:.1}K");
        }
    }

    #[test]
    fn debug_worst_case_species() {
        // Compare He/H at multiple altitudes against C reference
        let base = Nrlmsise00Input {
            day_of_year: 172,
            ut_seconds: 43200.0,
            altitude_km: 0.0,
            latitude_deg: 75.0,
            longitude_deg: 270.0,
            local_solar_time_hours: 12.0 + 270.0 / 15.0,
            f107_daily: 70.0,
            f107_avg: 70.0,
            ap_daily: 4.0,
            ap_array: [4.0; 7],
        };
        eprintln!("He/H density (lat=75, summer_solstice, solar_min):");
        let c_he = [
            3.181456e6, 7.021742e5, 3.983057e5, 2.398846e5, 1.470449e5, 5.761877e4, 1.554627e4,
        ];
        let c_h = [
            3.004045e6, 2.063095e5, 1.521022e5, 1.335147e5, 1.181220e5, 9.345534e4, 6.735497e4,
        ];
        for (i, &alt) in [120.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0]
            .iter()
            .enumerate()
        {
            let mut inp = base.clone();
            inp.altitude_km = alt;
            let (d, _, temp) = compute(&inp, &[1.0; 24]);
            let he_err = (d[0] - c_he[i]) / c_he[i] * 100.0;
            let h_err = (d[6] - c_h[i]) / c_h[i] * 100.0;
            eprintln!(
                "  {alt:7.1}km: He={:.6e} ({he_err:+.2}%)  H={:.6e} ({h_err:+.2}%)  T={temp:.2}",
                d[0], d[6]
            );
        }

        // Equatorial reference
        let base2 = Nrlmsise00Input {
            day_of_year: 80,
            ut_seconds: 43200.0,
            altitude_km: 0.0,
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            local_solar_time_hours: 12.0,
            f107_daily: 150.0,
            f107_avg: 150.0,
            ap_daily: 15.0,
            ap_array: [15.0; 7],
        };
        eprintln!("\nHe density (equatorial, vernal_equinox, solar_mod):");
        let c_he2 = [
            4.326096e7, 1.347404e7, 8.282984e6, 5.650467e6, 3.943564e6, 1.987947e6, 7.635356e5,
        ];
        for (i, &alt) in [120.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0]
            .iter()
            .enumerate()
        {
            let mut inp = base2.clone();
            inp.altitude_km = alt;
            let (d, _, temp) = compute(&inp, &[1.0; 24]);
            let he_err = (d[0] - c_he2[i]) / c_he2[i] * 100.0;
            eprintln!(
                "  {alt:7.1}km: He={:.6e} ({he_err:+.2}%)  T={temp:.2}",
                d[0]
            );
        }
    }

    #[test]
    fn composition_correction_dual_limits() {
        // High altitude: should return 1.0
        let v = composition_correction_dual(200.0, 0.5, 10.0, 50.0, 5.0);
        assert!(
            (v - 1.0).abs() < 1e-6,
            "composition_correction_dual high alt: {v}"
        );

        // Very low altitude (both e < -70): should return exp(r)
        let v = composition_correction_dual(-800.0, 0.5, 10.0, 50.0, 5.0);
        assert!(
            (v - 0.5_f64.exp()).abs() < 1e-10,
            "composition_correction_dual low alt: {v}"
        );

        // One e < -70, other normal: should be finite
        let v = composition_correction_dual(30.0, 0.5, 10.0, 50.0, 0.1);
        assert!(v.is_finite(), "composition_correction_dual mixed: {v}");
    }
}
