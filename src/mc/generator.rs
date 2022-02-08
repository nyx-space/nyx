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
extern crate rand_distr;

use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::{trajectory::InterpState, StateParameter};
use crate::{NyxError, Orbit};
use rand_distr::{Distribution, Normal};
use rand_pcg::Pcg64Mcg;

/// A state generator for Monte Carlo analyses.
pub struct Generator<S: InterpState, D: Distribution<f64> + Copy>
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
    pub dispersions: Vec<(StateParameter, D)>,
}

impl<S: InterpState, D: Distribution<f64> + Copy> Generator<S, D>
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
    pub fn new(
        template: S,
        params: &[StateParameter],
        dispersions: &[D],
    ) -> Result<Self, NyxError> {
        let mut me: Self = template.into();
        for (param, dispersion) in params.iter().zip(dispersions) {
            me.add_dispersion(*param, *dispersion)?;
        }
        Ok(me)
    }
}

impl<S: InterpState, D: Distribution<f64> + Copy> From<S> for Generator<S, D>
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

impl<S: InterpState> Generator<S, Normal<f64>>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
{
    /// Add a state dispersion from the provided 3-sigma value, zero mean
    pub fn add_from_3σ(
        &mut self,
        param: StateParameter,
        three_sigma: f64,
    ) -> Result<(), NyxError> {
        self.add_from_1σ(param, three_sigma / 3.0)
    }

    /// Add a state dispersion from the provided 1-sigma value, zero mean
    pub fn add_from_1σ(&mut self, param: StateParameter, std_dev: f64) -> Result<(), NyxError> {
        match self.template.value(&param) {
            Ok(_) => {
                self.dispersions
                    .push((param, Normal::new(0.0, std_dev).unwrap()));
                Ok(())
            }
            Err(e) => Err(e),
        }
    }
}

impl<S: InterpState, D: Distribution<f64> + Copy> Distribution<S> for Generator<S, D>
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
/// Generates any number of orbits from a template for any distribution D
pub type OrbitGenerator<D> = Generator<Orbit, D>;
/// Generates any number of spacecraft from a template for any distribution D
pub type SpacecraftGenerator<D> = Generator<Spacecraft, D>;
/// Generates any number of orbits from a template using a Normal distribution on the specified parameters
pub type GaussianOrbitGenerator = Generator<Orbit, Normal<f64>>;
/// Generates any number of spacecraft from a template using a Normal distribution on the specified parameters
pub type GaussianSpacecraftGenerator = Generator<Spacecraft, Normal<f64>>;

#[test]
fn generate_orbit() {
    use crate::cosmic::{Cosm, Orbit};
    use crate::time::Epoch;
    let cosm = Cosm::de438();

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    let mut orbit_generator: GaussianOrbitGenerator = state.into();
    orbit_generator
        .add_from_1σ(StateParameter::SMA, 1.0)
        .unwrap();

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

    // Allow for 10% error, which is huge
    assert!(
        cnt_too_far < 380,
        "Should have less than 33% of samples being more than 1 sigma away, got {}",
        cnt_too_far
    );
}
