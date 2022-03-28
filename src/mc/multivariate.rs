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

use std::iter::zip;

use super::rand_distr::{Distribution, Normal};
use super::DispersedState;
use crate::linalg::allocator::Allocator;
use crate::linalg::{
    Const, DefaultAllocator, DimMin, DimMinimum, DimName, DimSub, OMatrix, OVector,
};
use crate::md::StateParameter;
use crate::{NyxError, State};

/// A state generator for Monte Carlo analyses.
pub struct MultivariateNormal<S: State>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <<S::Size as DimMin<S::Size>>::Output as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output, S::Size>
        + Allocator<f64, <S::Size as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimSub<Const<1>>>::Output>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
    S::Size: DimMin<S::Size>,
    <S::Size as DimMin<S::Size>>::Output: DimSub<Const<1>>,
    S::Size: DimSub<Const<1>>,
{
    /// The template state
    pub template: S,
    /// The ordered vector of parameters to which the mean and covariance correspond to.
    pub params: Vec<StateParameter>,
    /// The mean of the multivariate normal distribution
    pub mean: OVector<f64, DimMinimum<S::Size, S::Size>>,
    /// The dot product \sqrt{\vec s} \cdot \vec v, where S is the singular values and V the V matrix from the SVD decomp of the covariance of multivariate normal distribution
    // pub sqrt_s_v: OVector<f64, DimMinimum<S::Size, S::Size>>,
    pub sqrt_s_v: OMatrix<f64, S::Size, DimMinimum<S::Size, S::Size>>,
    /// The standard normal distribution used to seed the multivariate normal distribution
    pub std_norm_distr: Normal<f64>,
}

impl<S: State> MultivariateNormal<S>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <<S::Size as DimMin<S::Size>>::Output as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output, S::Size>
        + Allocator<f64, <S::Size as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimSub<Const<1>>>::Output>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
    S::Size: DimMin<S::Size>,
    <S::Size as DimMin<S::Size>>::Output: DimSub<Const<1>>,
    S::Size: DimSub<Const<1>>,
{
    /// Creates a new Monte Carlos state generator from a mean and covariance which must be of the same size as the state vector
    /// The covariance must be positive semi definite. The algorithm is the one from numpy
    /// <https://github.com/numpy/numpy/blob/6c16f23c30fe490422959d30c2e22345211a2fe3/numpy/random/mtrand.pyx#L3979>
    pub fn new(
        template: S,
        params: Vec<StateParameter>,
        mean: OVector<f64, DimMinimum<S::Size, S::Size>>,
        cov: OMatrix<f64, S::Size, S::Size>,
    ) -> Result<Self, NyxError> {
        // Check that covariance is PSD by ensuring that all the eigenvalues are positive or nil
        match cov.eigenvalues() {
            None => return Err(NyxError::CovarianceMatrixNotPsd),
            Some(evals) => {
                for eigenval in &evals {
                    if *eigenval < 0.0 {
                        return Err(NyxError::CovarianceMatrixNotPsd);
                    }
                }
            }
        };

        let svd = cov.svd_unordered(false, true);
        if svd.v_t.is_none() {
            return Err(NyxError::CovarianceMatrixNotPsd);
        }

        let sqrt_s = svd.singular_values.map(|x| x.sqrt());
        let mut sqrt_s_v_t = svd.v_t.unwrap();

        for (i, mut col) in sqrt_s_v_t.column_iter_mut().enumerate() {
            col *= sqrt_s[i];
        }

        // OMatrix<f64, S::Size, <S::Size as DimMin<S::Size>>::Output, <DefaultAllocator as Allocator<f64, S::Size, <S::Size as DimMin<S::Size>>::Output>>::Buffer>

        Ok(Self {
            template,
            params,
            mean,
            sqrt_s_v: sqrt_s_v_t.transpose(),
            std_norm_distr: Normal::new(0.0, 1.0).unwrap(),
        })
    }

    /// Same as `new` but with a zero mean
    pub fn zero_mean(
        template: S,
        params: Vec<StateParameter>,
        cov: OMatrix<f64, S::Size, S::Size>,
    ) -> Result<Self, NyxError>
    where
        <S::Size as DimMin<S::Size>>::Output: DimName,
    {
        Self::new(
            template,
            params,
            OVector::<f64, DimMinimum<S::Size, S::Size>>::zeros(),
            cov,
        )
    }
}

impl<S: State> Distribution<DispersedState<S>> for MultivariateNormal<S>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <<S::Size as DimMin<S::Size>>::Output as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimMin<S::Size>>::Output>
        + Allocator<f64, <S::Size as DimMin<S::Size>>::Output, S::Size>
        + Allocator<f64, <S::Size as DimSub<Const<1>>>::Output>
        + Allocator<f64, S::Size, <S::Size as DimSub<Const<1>>>::Output>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
    S::Size: DimMin<S::Size>,
    <S::Size as DimMin<S::Size>>::Output: DimSub<Const<1>>,
    S::Size: DimSub<Const<1>>,
    <S::Size as DimMin<S::Size>>::Output: DimName,
{
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> DispersedState<S> {
        // Generate the vector representing the state
        let x_rng = OVector::<f64, S::Size>::from_fn(|_, _| self.std_norm_distr.sample(rng));
        let x = self.sqrt_s_v.transpose() * x_rng + &self.mean;
        let mut state = self.template;

        let mut actual_dispersions = Vec::new();
        for (delta, param) in zip(&x, &self.params) {
            actual_dispersions.push((*param, *delta));
            // We know this state can return something for this param
            let cur_value = state.value(&param).unwrap();
            state.set_value(&param, cur_value + delta).unwrap();
        }

        DispersedState {
            state,
            actual_dispersions,
        }
    }
}

#[test]
fn test_multivariate_state() {
    use crate::cosmic::{Cosm, Orbit};
    use crate::linalg::{Matrix6, Vector6};
    use crate::time::Epoch;
    use rand_pcg::Pcg64Mcg;

    let cosm = Cosm::de438();

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    let mean = Vector6::zeros();
    let std_dev = Vector6::new(10.0, 10.0, 10.0, 0.2, 0.2, 0.2);
    let cov = Matrix6::from_diagonal(&std_dev);

    let orbit_generator = state.disperse(mean, cov).unwrap();

    // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
    // Create a reproducible fast seed
    let seed = 0;
    let rng = Pcg64Mcg::new(seed);

    let cnt_too_far: u16 = orbit_generator
        .sample_iter(rng)
        .take(1000)
        .map(|dispersed_state| {
            let mut cnt = 0;
            for idx in 0..6 {
                let val_std_dev = std_dev[idx];
                let cur_val = dispersed_state.state.as_vector().unwrap()[idx];
                let nom_val = state.as_vector().unwrap()[idx];
                if (cur_val - nom_val).abs() > val_std_dev {
                    cnt += 1;
                }
            }
            cnt
        })
        .sum::<u16>();

    // We specified a seed so we know exactly what to expect
    assert_eq!(
        cnt_too_far / 6,
        329,
        "Should have less than 33% of samples being more than 1 sigma away, got {}",
        cnt_too_far
    );
}
