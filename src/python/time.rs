/// nyx_space.time
use pyo3::prelude::*;
use pyo3::types::PyType;
use pyo3::exceptions::PyException;

use crate::time::Epoch as EpochRs;

#[pyclass]
pub struct Epoch {
    pub inner: EpochRs,
}

// Required as hifitime is now a separate crate
#[derive(Debug)]
pub struct TimeError(crate::time::Errors);

impl std::convert::From<crate::time::Errors> for TimeError
{
    fn from(other: crate::time::Errors) -> Self {
        Self(other)
    }
}

impl std::convert::From<TimeError> for PyErr 
{
    fn from(err: TimeError) -> Self {
        PyException::new_err(err.0.to_string())
    }
}


#[pymethods]
impl Epoch {
    /// Initialize an Epoch from the provided TAI seconds since 1900 January 01 at midnight
    #[classmethod]
    pub fn from_tai_seconds(_cls: &PyType, seconds: f64) -> Self {
        Self {
            inner: EpochRs::from_tai_seconds(seconds)
        }
    }

    /// Initialize an Epoch from the provided TAI days since 1900 January 01 at midnight
    #[classmethod]
    pub fn from_tai_days(_cls: &PyType, days: f64) -> Self {
        Self {
            inner: EpochRs::from_tai_days(days)
        }
    }

    #[classmethod]
    pub fn from_mjd_tai(_cls: &PyType, days: f64) -> Self {
        Self { inner: EpochRs::from_mjd_tai(days) }
    }

    #[classmethod]
    pub fn from_jde_tai(_cls: &PyType, days: f64) -> Self {
        Self {
            inner: EpochRs::from_jde_tai(days)
        }
    }

    /// Initialize an Epoch from the provided TT seconds (approximated to 32.184s delta from TAI)
    #[classmethod]
    pub fn from_tt_seconds(_cls: &PyType, seconds: f64) -> Self {
        Self {
            inner: EpochRs::from_tt_seconds(seconds)
        }
    }

    #[classmethod]
    pub fn from_et_seconds(_cls: &PyType, seconds: f64) -> Epoch {
        Self {
            inner: EpochRs::from_et_seconds(seconds)
        }
    }

    /// Initialize from Dynamic Barycentric Time (TDB) (same as SPICE ephemeris time) whose epoch is 2000 JAN 01 noon TAI
    #[classmethod]
    pub fn from_tdb_seconds(_cls: &PyType, seconds: f64) -> Epoch {
        Self {
            inner: EpochRs::from_tdb_seconds(seconds)
        }
    }

    // /// Initialize from Dynamic Barycentric Time (TDB) (same as SPICE ephemeris time) whose epoch is 2000 JAN 01 noon TAI
    // fn from_tdb_seconds_d(duration: Duration) -> Epoch {
    //     Self {
    //         inner: EpochRs::from_tdb_seconds_d(duration)
    //     }
    // }

    #[classmethod]
    pub fn from_jde_et(_cls: &PyType, days: f64) -> Self {
        Self {
            inner: EpochRs::from_jde_et(days)
        }
    }

    /// Initialize from Dynamic Barycentric Time (TDB) (same as SPICE ephemeris time) in JD days
    #[classmethod]
    pub fn from_jde_tdb(_cls: &PyType, days: f64) -> Self {
        Self {
            inner: EpochRs::from_jde_tdb(days)
        }
    }

    // /// Builds a new Epoch from the hi and lo two-float values
    // #[classmethod]
    // pub fn try_from_hi_lo(_cls: &PyType, hi: f64, lo: f64) -> PyResult<Self> {
    //     Ok(Self {
    //         inner: EpochRs::try_from_hi_lo(hi, lo)?
    //     })
    // }

    /// Attempts to build an Epoch from the provided Gregorian date and time in TAI.
    #[classmethod]
    pub fn maybe_from_gregorian_tai(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: u32,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: EpochRs::maybe_from_gregorian_tai(year, month, day, hour, minute, second, nanos).map_err(TimeError::from)?
        })
    }

    // /// Attempts to build an Epoch from the provided Gregorian date and time in the provided time system.
    // #[allow(clippy::too_many_arguments)]
    // pub fn maybe_from_gregorian(
    //     year: i32,
    //     month: u8,
    //     day: u8,
    //     hour: u8,
    //     minute: u8,
    //     second: u8,
    //     nanos: u32,
    //     ts: TimeSystem,
    // ) -> PyResult<Self> {
        
    // }

    /// Builds an Epoch from the provided Gregorian date and time in TAI. If invalid date is provided, this function will panic.
    /// Use maybe_from_gregorian_tai if unsure.
    #[classmethod]
    pub fn from_gregorian_tai(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: u32,
    ) -> Self {
        Self {
            inner: EpochRs::from_gregorian_tai(year, month, day, hour, minute, second, nanos)
        }
    }

    #[classmethod]
    pub fn from_gregorian_tai_at_midnight(_cls: &PyType, year: i32, month: u8, day: u8) -> Self {
        Self {
            inner: EpochRs::from_gregorian_tai_at_midnight(year, month, day)
        }
    }

    #[classmethod]
    pub fn from_gregorian_tai_at_noon(_cls: &PyType, year: i32, month: u8, day: u8) -> Self {
        Self {
            inner: EpochRs::from_gregorian_tai_at_noon(year, month, day)
        }
    }

    #[classmethod]
    pub fn from_gregorian_tai_hms(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
    ) -> Self {
        Self {
            inner: EpochRs::from_gregorian_tai_hms(year, month, day, hour, minute, second)
        }
    }

    /// Attempts to build an Epoch from the provided Gregorian date and time in UTC.
    #[classmethod]
    pub fn maybe_from_gregorian_utc(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: u32,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: EpochRs::maybe_from_gregorian_utc(year, month, day, hour, minute, second, nanos).map_err(TimeError::from)?
        })
    }

    /// Builds an Epoch from the provided Gregorian date and time in TAI. If invalid date is provided, this function will panic.
    /// Use maybe_from_gregorian_tai if unsure.
    #[classmethod]
    pub fn from_gregorian_utc(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: u32,
    ) -> Self {
        Self {
            inner: EpochRs::from_gregorian_utc(year, month, day, hour, minute, second, nanos)
        }
    }

    #[classmethod]
    pub fn from_gregorian_utc_at_midnight(_cls: &PyType, year: i32, month: u8, day: u8) -> Self {
        Self {
            inner: EpochRs::from_gregorian_utc_at_midnight(year, month, day)
        }
    }

    #[classmethod]
    pub fn from_gregorian_utc_at_noon(_cls: &PyType, year: i32, month: u8, day: u8) -> Self {
        Self {
            inner: EpochRs::from_gregorian_utc_at_noon(year, month, day)
        }
    }

    #[classmethod]
    pub fn from_gregorian_utc_hms(
        _cls: &PyType,
        year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
    ) -> Self {
        Self {
            inner: EpochRs::from_gregorian_utc_hms(year, month, day, hour, minute, second)
        }
    }

    pub fn as_tai_seconds(&self) -> f64 {
        self.inner.as_tai_seconds()
    }

    // /// Returns this time in a Duration past J1900 counted in TAI
    // pub fn as_tai_duration(&self) -> Duration {
    //     self.inner.as_tai_seconds()
    // }

    // /// Returns the epoch as a floating point value in the provided unit
    // pub fn as_tai(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_tai()
    // }

    pub fn as_tai_days(&self) -> f64 {
        self.inner.as_tai_days()
    }

    /// Returns the number of UTC seconds since the TAI epoch
    pub fn as_utc_seconds(&self) -> f64 {
        self.inner.as_utc_seconds()
    }

    // /// Returns this time in a Duration past J1900 counted in UTC
    // fn as_utc_duration(&self) -> Duration {
    //     self.inner.as_utc_seconds()
    // }

    // /// Returns the number of UTC seconds since the TAI epoch
    // pub fn as_utc(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_utc()
    // }

    /// Returns the number of UTC days since the TAI epoch
    pub fn as_utc_days(&self) -> f64 {
        self.inner.as_utc_days()
    }

    /// `as_mjd_days` creates an Epoch from the provided Modified Julian Date in days as explained
    /// [here](http://tycho.usno.navy.mil/mjd.html). MJD epoch is Modified Julian Day at 17 November 1858 at midnight.
    pub fn as_mjd_tai_days(&self) -> f64 {
        self.inner.as_mjd_tai_days()
    }

    /// Returns the Modified Julian Date in seconds TAI.
    pub fn as_mjd_tai_seconds(&self) -> f64 {
        self.inner.as_mjd_tai_seconds()
    }

    // pub fn as_mjd_tai(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_mjd_tai(unit)
    // }

    /// Returns the Modified Julian Date in days UTC.
    pub fn as_mjd_utc_days(&self) -> f64 {
        self.inner.as_mjd_utc_days()
    }

    // /// Returns the Modified Julian Date in the provided unit in UTC.
    // pub fn as_mjd_utc(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_mjd_utc_days()
    // }

    /// Returns the Modified Julian Date in seconds UTC.
    pub fn as_mjd_utc_seconds(&self) -> f64 {
        self.inner.as_mjd_utc_seconds()
    }

    /// Returns the Julian days from epoch 01 Jan -4713, 12:00 (noon)
    /// as explained in "Fundamentals of astrodynamics and applications", Vallado et al.
    /// 4th edition, page 182, and on [Wikipedia](https://en.wikipedia.org/wiki/Julian_day).
    pub fn as_jde_tai_days(&self) -> f64 {
        self.inner.as_jde_tai_days()
    }

    // pub fn as_jde_tai(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_jde_tai()
    // }

    // pub fn as_jde_tai_duration(&self) -> Duration {
    //     self.inner.as_jde_tai_duration()
    // }

    /// Returns the Julian seconds in TAI.
    pub fn as_jde_tai_seconds(&self) -> f64 {
        self.inner.as_jde_tai_seconds()
    }

    /// Returns the Julian days in UTC.
    pub fn as_jde_utc_days(&self) -> f64 {
        self.inner.as_jde_utc_days()
    }

    // pub fn as_jde_utc_duration(&self) -> Duration {
    //     self.as_utc_duration() + TimeUnit::Day * (J1900_OFFSET + MJD_OFFSET)
    // }

    /// Returns the Julian seconds in UTC.
    pub fn as_jde_utc_seconds(&self) -> f64 {
        self.inner.as_jde_utc_seconds()
    }

    /// Returns seconds past TAI epoch in Terrestrial Time (TT) (previously called Terrestrial Dynamical Time (TDT))
    pub fn as_tt_seconds(&self) -> f64 {
        self.inner.as_tt_seconds()
    }

    // pub fn as_tt_duration(&self) -> Duration {
    //     self.0 + TimeUnit::Second * TT_OFFSET_S
    // }

    /// Returns days past TAI epoch in Terrestrial Time (TT) (previously called Terrestrial Dynamical Time (TDT))
    pub fn as_tt_days(&self) -> f64 {
        self.inner.as_tt_days()
    }

    /// Returns the centuries pased J2000 TT
    pub fn as_tt_centuries_j2k(&self) -> f64 {
        self.inner.as_tt_centuries_j2k()
    }

    // /// Returns the duration past J2000 TT
    // pub fn as_tt_since_j2k(&self) -> Duration {
    //     self.as_tt_duration() - TimeUnit::Second * ET_EPOCH_S
    // }

    /// Returns days past Julian epoch in Terrestrial Time (TT) (previously called Terrestrial Dynamical Time (TDT))
    pub fn as_jde_tt_days(&self) -> f64 {
        self.inner.as_jde_tt_days()
    }

    // pub fn as_jde_tt_duration(&self) -> Duration {
    //     self.as_tt_duration() + TimeUnit::Day * (J1900_OFFSET + MJD_OFFSET)
    // }

    /// Returns days past Modified Julian epoch in Terrestrial Time (TT) (previously called Terrestrial Dynamical Time (TDT))
    pub fn as_mjd_tt_days(&self) -> f64 {
        self.inner.as_mjd_tt_days()
    }

    // pub fn as_mjd_tt_duration(&self) -> Duration {
    //     self.as_tt_duration() + TimeUnit::Day * J1900_OFFSET
    // }

    /// Returns seconds past GPS Time Epoch, defined as UTC midnight of January 5th to 6th 1980 (cf. https://gssc.esa.int/navipedia/index.php/Time_References_in_GNSS#GPS_Time_.28GPST.29).
    pub fn as_gpst_seconds(&self) -> f64 {
        self.inner.as_gpst_seconds()
    }

    // pub fn as_gpst_duration(&self) -> Duration {
    //     self.as_tai_duration() - TimeUnit::Second * 19.0
    // }

    /// Returns days past GPS Time Epoch, defined as UTC midnight of January 5th to 6th 1980 (cf. https://gssc.esa.int/navipedia/index.php/Time_References_in_GNSS#GPS_Time_.28GPST.29).
    pub fn as_gpst_days(&self) -> f64 {
        self.inner.as_gpst_days()
    }

    /// Returns the Ephemeris Time seconds past epoch
    pub fn as_et_seconds(&self) -> f64 {
        self.inner.as_et_seconds()
    }

    // pub fn as_et_duration(&self) -> Duration {
    //     self.as_tai_duration() - TimeUnit::Second * (ET_EPOCH_S - ET_OFFSET_S)
    // }

    /// Returns the Dynamic Barycentric Time (TDB) (higher fidelity SPICE ephemeris time) whose epoch is 2000 JAN 01 noon TAI (cf. https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TDT_-_TDB.2C_TCB)
    pub fn as_tdb_seconds(&self) -> f64 {
        self.inner.as_tdb_seconds()
    }

    // pub fn as_tdb_duration(&self) -> Duration {
    //     let inner = self.inner_g_rad();

    //     self.as_tt_duration() - TimeUnit::Second * (ET_EPOCH_S + (0.001_658 * inner.sin()))
    // }

    /// Returns the Ephemeris Time JDE past epoch
    pub fn as_jde_et_days(&self) -> f64 {
        self.inner.as_jde_et_days()
    }

    // pub fn as_jde_et_duration(&self) -> Duration {
    //     self.as_jde_tt_duration() + TimeUnit::Second * 0.000_935
    // }

    // pub fn as_jde_et(self, unit: TimeUnit) -> f64 {
    //     self.inner.as_jde_et()
    // }

    // pub fn as_jde_tdb_duration(&self) -> Duration {
    //     self.as_jde_tdb_days() * TimeUnit::Day
    // }

    /// Returns the Dynamic Barycentric Time (TDB) (higher fidelity SPICE ephemeris time) whose epoch is 2000 JAN 01 noon TAI (cf. https://gssc.esa.int/navipedia/index.php/Transformations_between_Time_Systems#TDT_-_TDB.2C_TCB)
    pub fn as_jde_tdb_days(&self) -> f64 {
        self.inner.as_jde_tdb_days()
    }

    /// Returns the number of days since Dynamic Barycentric Time (TDB) J2000 (used for Archinal et al. rotations)
    pub fn as_tdb_days_since_j2000(&self) -> f64 {
        self.inner.as_tdb_days_since_j2000()
    }

    /// Returns the number of centuries since Dynamic Barycentric Time (TDB) J2000 (used for Archinal et al. rotations)
    pub fn as_tdb_centuries_since_j2000(&self) -> f64 {
        self.inner.as_tdb_centuries_since_j2000()
    }

    /// Converts an ISO8601 Datetime representation without timezone offset to an Epoch.
    /// If no time system is specified, than UTC is assumed.
    /// The `T` which separates the date from the time can be replaced with a single whitespace character (`\W`).
    /// The offset is also optional, cf. the examples below.
    ///
    /// # Example
    /// ```
    /// use hifitime::Epoch;
    /// let dt = Epoch::from_gregorian_utc(2017, 1, 14, 0, 31, 55, 0);
    /// assert_eq!(
    ///     dt,
    ///     Epoch::from_gregorian_str("2017-01-14T00:31:55 UTC").unwrap()
    /// );
    /// assert_eq!(
    ///     dt,
    ///     Epoch::from_gregorian_str("2017-01-14T00:31:55.0000 UTC").unwrap()
    /// );
    /// assert_eq!(
    ///     dt,
    ///     Epoch::from_gregorian_str("2017-01-14T00:31:55").unwrap()
    /// );
    /// assert_eq!(
    ///     dt,
    ///     Epoch::from_gregorian_str("2017-01-14 00:31:55").unwrap()
    /// );
    /// // Regression test for #90
    /// assert_eq!(
    ///     Epoch::from_gregorian_utc(2017, 1, 14, 0, 31, 55, 811000000),
    ///     Epoch::from_gregorian_str("2017-01-14 00:31:55.811 UTC").unwrap()
    /// );
    /// assert_eq!(
    ///     Epoch::from_gregorian_utc(2017, 1, 14, 0, 31, 55, 811200000),
    ///     Epoch::from_gregorian_str("2017-01-14 00:31:55.8112 UTC").unwrap()
    /// );
    /// ```
    #[classmethod]
    pub fn from_gregorian_str(_cls: &PyType, s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: EpochRs::from_gregorian_str(s).map_err(TimeError::from)?
        })
    }

    /// Converts the Epoch to the Gregorian UTC equivalent as (year, month, day, hour, minute, second).
    /// WARNING: Nanoseconds are lost in this conversion!
    ///
    /// # Example
    /// ```
    /// use hifitime::Epoch;
    /// let dt_str = "2017-01-14T00:31:55 UTC";
    /// let dt = Epoch::from_gregorian_str(dt_str).unwrap();
    /// let (y, m, d, h, min, s, _) = dt.as_gregorian_utc();
    /// assert_eq!(y, 2017);
    /// assert_eq!(m, 1);
    /// assert_eq!(d, 14);
    /// assert_eq!(h, 0);
    /// assert_eq!(min, 31);
    /// assert_eq!(s, 55);
    /// assert_eq!(dt_str, dt.as_gregorian_utc_str().to_owned());
    /// ```
    pub fn as_gregorian_utc(&self) -> (i32, u8, u8, u8, u8, u8, u32) {
        self.inner.as_gregorian_utc()
    }

    /// Converts the Epoch to UTC Gregorian in the ISO8601 format.
    pub fn as_gregorian_utc_str(&self) -> String {
        self.inner.as_gregorian_utc_str()
    }

    /// Converts the Epoch to the Gregorian TAI equivalent as (year, month, day, hour, minute, second).
    /// WARNING: Nanoseconds are lost in this conversion!
    ///
    /// # Example
    /// ```
    /// use hifitime::Epoch;
    /// let dt = Epoch::from_gregorian_tai_at_midnight(1972, 1, 1);
    /// let (y, m, d, h, min, s, _) = dt.as_gregorian_tai();
    /// assert_eq!(y, 1972);
    /// assert_eq!(m, 1);
    /// assert_eq!(d, 1);
    /// assert_eq!(h, 0);
    /// assert_eq!(min, 0);
    /// assert_eq!(s, 0);
    /// ```
    pub fn as_gregorian_tai(&self) -> (i32, u8, u8, u8, u8, u8, u32) {
        self.inner.as_gregorian_tai()
    }

    /// Converts the Epoch to TAI Gregorian in the ISO8601 format with " TAI" appended to the string
    pub fn as_gregorian_tai_str(&self) -> String {
        self.inner.as_gregorian_tai_str()
    }

    // /// Converts the Epoch to Gregorian in the provided time system and in the ISO8601 format with the time system appended to the string
    // pub fn as_gregorian_str(self, ts: TimeSystem) -> String {
        
    // }
}

#[pymodule]
pub fn time(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Epoch>()?;
    Ok(())
}
