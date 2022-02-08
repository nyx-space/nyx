/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::rand_distr::{Distribution, Normal};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::StateParameter;
use crate::{NyxError, State};

/// A state generator for Monte Carlo analyses.
pub struct Generator<S: State, Distr: Distribution<f64> + Copy>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    /// The template state
    pub template: S,
    /// The list of dispersions to be added to the template state
    /// Note that we can't use a HashMap here because StateParameter has a SlantAngle option comprised of f64s, and those neither have Hash nor Eq
    pub dispersions: Vec<(StateParameter, Distr)>,
}

impl<S: State, D: Distribution<f64> + Copy> Generator<S, D>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    /// Add a parameter dispersion to this Monte Carlo state generator.
    pub fn add_dispersion(&mut self, param: StateParameter, dispersion: D) -> Result<(), NyxError> {
        // Try to set that parameter, and report an error on initialization if it fails
        match self.template.clone().set_value(&param, 0.0) {
            Ok(_) => {
                self.dispersions.push((param, dispersion));
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    /// Create a new Monte Carlo state generator given a template state, the parameters to disperse, and their respective dispersion probability density functions.
    pub fn from_dispersions(
        template: S,
        dispersions: &[(StateParameter, D)],
    ) -> Result<Self, NyxError> {
        let mut me: Self = template.into();
        for (param, dispersion) in dispersions {
            me.add_dispersion(*param, *dispersion)?;
        }
        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameters to disperse, and their respective dispersion probability density functions.
    pub fn from_dispersion(
        template: S,
        param: StateParameter,
        dispersion: D,
    ) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        me.add_dispersion(param, dispersion)?;

        Ok(me)
    }
}

impl<S: State, D: Distribution<f64> + Copy> From<S> for Generator<S, D>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    fn from(template: S) -> Self {
        Self {
            template,
            dispersions: Vec::new(),
        }
    }
}

impl<S: State> Generator<S, Normal<f64>>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    /// Add a state dispersion from the provided 3-sigma value, zero mean
    pub fn add_3σ(&mut self, param: StateParameter, three_sigma: f64) -> Result<(), NyxError> {
        self.add_1σ(param, three_sigma / 3.0)
    }

    /// Add a state dispersion from the provided 1-sigma value, zero mean
    pub fn add_1σ(&mut self, param: StateParameter, std_dev: f64) -> Result<(), NyxError> {
        match self.template.value(&param) {
            Ok(_) => {
                self.dispersions
                    .push((param, Normal::new(0.0, std_dev).unwrap()));
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    /// Create a new Monte Carlo state generator given a template state, the parameters to disperse, and their respective 3-σ standard deviations, zero mean.
    pub fn from_3σs(
        template: S,
        three_sigmas: &[(StateParameter, f64)],
    ) -> Result<Self, NyxError> {
        let mut me: Self = template.into();
        for (param, three_sigma) in three_sigmas {
            me.add_3σ(*param, *three_sigma)?;
        }
        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameter to disperse, and its 3-σ standard deviation, zero mean.
    pub fn from_3σ(
        template: S,
        param: StateParameter,
        three_sigma: f64,
    ) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        me.add_3σ(param, three_sigma)?;

        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameter to disperse, and its 1-σ standard deviation in percentage of the template's value, zero mean.
    pub fn from_3σ_prct(template: S, param: StateParameter, prct: f64) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        me.add_3σ(param, template.value(&param)? * prct)?;

        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameters to disperse, and their respective 1-σ standard deviations, zero mean.
    pub fn from_1σs(template: S, std_devs: &[(StateParameter, f64)]) -> Result<Self, NyxError> {
        let mut me: Self = template.into();
        for (param, three_sigma) in std_devs {
            me.add_1σ(*param, *three_sigma)?;
        }
        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameter to disperse, and its 1-σ standard deviation in percentage of the template's value, zero mean.
    pub fn from_1σs_prct(template: S, prcts: &[(StateParameter, f64)]) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        for (param, prct) in prcts {
            me.add_1σ(*param, template.value(&param)? * prct)?;
        }

        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameter to disperse, and its 1-σ standard deviation, zero mean.
    pub fn from_1σ(template: S, param: StateParameter, std_dev: f64) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        me.add_1σ(param, std_dev)?;

        Ok(me)
    }

    /// Create a new Monte Carlo state generator given a template state, the parameter to disperse, and its 1-σ standard deviation in percentage of the template's value, zero mean.
    pub fn from_1σ_prct(template: S, param: StateParameter, prct: f64) -> Result<Self, NyxError> {
        let mut me: Self = template.into();

        me.add_1σ(param, template.value(&param)? * prct)?;

        Ok(me)
    }
}

impl<S: State, D: Distribution<f64> + Copy> Distribution<S> for Generator<S, D>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> S {
        let mut me = self.template;
        for (param, dispersion) in &self.dispersions {
            // We know this state can return something for this param
            let cur_value = me.value(param).unwrap();
            // Apply the dispersion
            let delta = dispersion.sample(rng);
            me.set_value(param, cur_value + delta).unwrap();
        }
        me
    }
}

/// Generates a state generator with a Normal distribution
pub type GaussianGenerator<S> = Generator<S, Normal<f64>>;

#[test]
fn generate_orbit() {
    use crate::cosmic::{Cosm, Orbit};
    use crate::time::Epoch;
    use rand_pcg::Pcg64Mcg;

    let cosm = Cosm::de438();

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    let orbit_generator = GaussianGenerator::from_1σ(state, StateParameter::SMA, 1.0).unwrap();

    // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
    // Create a reproducible fast seed
    let seed = 0;
    let rng = Pcg64Mcg::new(seed);

    let cnt_too_far: u16 = orbit_generator
        .sample_iter(rng)
        .take(1000)
        .map(|state| {
            if (8_191.93 - state.sma()).abs() > 1.0 {
                1
            } else {
                0
            }
        })
        .sum::<u16>();

    // We specified a seed so we know exactly what to expect
    assert_eq!(
        cnt_too_far, 308,
        "Should have less than 33% of samples being more than 1 sigma away, got {}",
        cnt_too_far
    );
}

#[test]
fn generate_spacecraft() {
    use crate::cosmic::{Cosm, Orbit, Spacecraft, State};
    use crate::dynamics::guidance::Thruster;
    use crate::time::Epoch;
    use rand_pcg::Pcg64Mcg;

    let cosm = Cosm::de438();

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let orbit = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    let nominal_thrust = 50.0;
    let nominal_isp = 300.0;

    let sc = Spacecraft::from_thruster(
        orbit,
        1_000.0,
        750.0,
        Thruster {
            isp_s: 300.0,
            thrust_N: 50.0,
        },
        crate::io::odp::GuidanceMode::Inhibit,
    );

    let sc_generator = GaussianGenerator::from_1σs_prct(
        sc,
        &[(StateParameter::Thrust, 0.05), (StateParameter::Isp, 0.01)],
    )
    .unwrap();

    // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
    // Create a reproducible fast seed
    let seed = 0;
    let rng = Pcg64Mcg::new(seed);

    let cnt_too_far: u16 = sc_generator
        .sample_iter(rng)
        .take(1000)
        .map(|state| {
            // Check out of bounds
            let thrust_oob = (nominal_thrust - state.value(&StateParameter::Thrust).unwrap()).abs()
                / nominal_thrust
                > 0.05;
            let isp_oob = (nominal_isp - state.value(&StateParameter::Isp).unwrap()).abs()
                / nominal_isp
                > 0.01;
            if thrust_oob || isp_oob {
                1
            } else {
                0
            }
        })
        .sum::<u16>();

    // We specified a seed, so we can specify exactly the number we're expecting
    assert_eq!(
        cnt_too_far, 511,
        "Should have less than 33% of samples two two draws being more than 1 sigma away, got {}",
        cnt_too_far
    );
}
