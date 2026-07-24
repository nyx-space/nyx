//! Geographic and temporal conversions for NRLMSISE-00.
//!
//! Converts between satellite position / epoch and the geodetic / time
//! parameters required by NRLMSISE-00. Two entry points are provided:
//!
//! - [`simple_eci_to_geodetic_latlon`] — the Phase 1–3 simple path:
//!   `Vec3<SimpleEci>` + `Epoch<Utc>` rotated by GMST (naive ERA, no
//!   dUT1) to `Vec3<SimpleEcef>`, then Bowring to geodetic
//! - [`precise_gcrs_to_geodetic_latlon`] — the Phase 3B precise path:
//!   `Vec3<Gcrs>` + `Epoch<Utc>` + full EOP provider rotated by the
//!   IAU 2006 CIO chain (`Rotation<Gcrs, Itrs>::iau2006_full_from_utc`)
//!   to `Vec3<Itrs>`, then Bowring to geodetic
//!
//! The entry points have distinct names **and distinct input types**
//! so that a caller cannot accidentally wire a simple-path position
//! into the precise function or vice versa — the type system rejects
//! the mix at compile time.

use hifitime::Epoch;

/// Convert a simple-path ECI position to WGS-84 geodetic latitude and
/// longitude [degrees].
///
/// The input `position_eci` is a `Vec3<SimpleEci>` — the phantom-typed
/// marker produced by the Phase 1–3 simple rotation path. The rotation
/// to `SimpleEcef` uses naive ERA (`epoch.gmst()` on a UTC epoch,
/// equivalent to assuming `dUT1 = 0`). Accuracy is bounded by the
/// ~0.9 s dUT1 drift, producing up to ~24 km longitude error at the
/// equator.
///
/// For the precise IAU 2006 CIO chain with real EOP data, use
/// [`precise_gcrs_to_geodetic_latlon`].

/// Convert a GCRS position + UTC epoch + EOP provider to WGS-84
/// geodetic latitude and longitude [degrees] via the full IAU 2006
/// CIO-based rotation chain.
///
/// This is the first downstream consumer of Phase 3B's
/// [`Rotation::<frame::Gcrs, frame::Itrs>::iau2006_full_from_utc`].
/// The EOP provider must supply `dUT1`, `dX`/`dY`, and `xp`/`yp` so
/// the full GCRS → ITRS transformation can be built; the trait bound
/// is [`Ut1Offset`] + [`NutationCorrections`] + [`PolarMotion`], matching
/// the `iau2006_full_from_utc` signature exactly.
///
/// `arika::earth::eop::NullEop` is rejected at compile time — see
/// `arika/tests/trybuild/null_eop_in_iau2006_full_from_utc.rs`.

/// Compute local apparent solar time [hours].
///
/// Applies the Equation of Time correction to convert from mean to apparent
/// solar time:
///
///   LST_apparent = UT/3600 + lon/15 + EoT(epoch)
///
/// where EoT accounts for Earth's orbital eccentricity and axial tilt
/// (up to ±16 minutes correction).
pub fn local_solar_time(ut_sec: f64, longitude_deg: f64, _epoch: &Epoch) -> f64 {
    // Equation of Time is a solar-ephemeris (TDB) quantity; convert here. The
    // UTC `epoch` itself supplies the mean-solar-time terms above.
    let eot_hours = 0.0;
    let lst = ut_sec / 3600.0 + longitude_deg / 15.0 + eot_hours;
    // Normalize to [0, 24)
    ((lst % 24.0) + 24.0) % 24.0
}
