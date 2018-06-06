extern crate hifitime;
extern crate nalgebra as na;

use self::rand::distributions::{IndependentSample, Normal};
use self::rand::thread_rng;
use super::Measurement;
use celestia::{CoordinateFrame, State, ECEF, ECI};
use hifitime::instant::Instant;
use na::{Matrix2x6, U2, U6, Vector2};

/// GroundStation defines a Two Way ranging equipment.
#[derive(Debug, Clone, Copy)]
pub struct GroundStation {
    pub name: String,
    pub elevation_mask: f64,
    range_noise: Normal,
    range_rate_noise: Normal,
}

impl GroundStation {
    /// Initializes a new Two Way ranging equipment from the noise values.
    pub fn from_noise_values(name: String, elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation {
            name,
            elevation_mask,
            range_noise: Normal::new(1.0, range_noise),
            range_rate_noise: Normal::new(1.0, range_rate_noise),
        }
    }

    /// Perform a measurement from the transmitter (tx) to the receiver (rx).
    pub fn measure(self, tx: State, rx: State, dt: Instant) -> GroundMeasurement {
        GroundMeasurement::new(
            tx,
            rx,
            self.elevation_mask,
            Vector2::new(
                self.range_noise.ind_sample(&mut thread_rng()),
                self.range_rate_noise.ind_sample(&mut thread_rng()),
            ),
        )
    }
}

/// Implements the Range and Range Rate measurement.
#[derive(Debug, Clone, Copy)]
pub struct GroundMeasurement {
    visible: bool,
    obs: Vector2<f64>,
    h_tilde: Matrix2x6<f64>,
    dt: Instant,
}

impl GroundMeasurement {
    pub fn range(&self) -> f64 {
        self.obs[(0, 0)]
    }
    pub fn range_rate(&self) -> f64 {
        self.obs[(1, 0)]
    }
}

impl Measurement for GroundMeasurement {
    type StateSize = U6;
    type MeasurementSize = U2;

    pub fn new(tx: State<ECI>, rx: State<ECI>, elevation_mask: f64, noise: Vector2<f64>) -> GroundMeasurement {
        // Let's start by computing the range and range rate
        // Ensure both are in the ECEF frame
        let tx_ecef = tx.in_frame(ECEF {});
        let rx_ecef = rx.in_frame(ECEF {});
        let rho_ecef = rx_ecef - tx_ecef;
        // Convert to SEZ frame to compute the elevation.

        /*
        ρECEF = make([]float64, 3)
	for i := 0; i < 3; i++ {
		ρECEF[i] = rECEF[i] - s.R[i]
	}
	ρ = Norm(ρECEF)
	rSEZ := MxV33(R3(s.Longθ), ρECEF)
	rSEZ = MxV33(R2(math.Pi/2-s.LatΦ), rSEZ)
	el = math.Asin(rSEZ[2]/ρ) * r2d
	az = (2*math.Pi + math.Atan2(rSEZ[1], -rSEZ[0])) * r2d
return*/

        /*
        // The station vectors are in ECEF, so let's convert the state to ECEF.
	rECEF := ECI2ECEF(state.Orbit.R(), θgst)
	vECEF := ECI2ECEF(state.Orbit.V(), θgst)
	// Compute visibility for each station.
	ρECEF, ρ, el, _ := s.RangeElAz(rECEF)
	vDiffECEF := make([]float64, 3)
	for i := 0; i < 3; i++ {
		vDiffECEF[i] = (vECEF[i] - s.V[i]) / ρ
	}
	ρDot := mat64.Dot(mat64.NewVector(3, ρECEF), mat64.NewVector(3, vDiffECEF))
	ρNoisy := ρ + s.RangeNoise.Rand(nil)[0]
ρDotNoisy := ρDot + s.RangeRateNoise.Rand(nil)[0]
*/
    }

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    pub fn observation(&self) -> &Vector2<f64> {
        &self.obs
    }

    pub fn sensitivity(&self) -> &Matrix2x6<f64> {
        &self.h_tilde
    }
}
