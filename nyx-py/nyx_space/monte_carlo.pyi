from __future__ import annotations
import typing

@typing.final
class MvnSpacecraft:
    """A multivariate spacecraft state generator for Monte Carlo analyses. Ensures that the covariance is properly applied on all provided state variables.

# Algorithm

The `MvnSpacecraft` allows sampling from a multivariate normal distribution defined by a template state (mean) and a set of state dispersions (uncertainties).
The core difficulty is that dispersions are often defined in non-Cartesian spaces (like Keplerian orbital elements or B-Plane parameters), while the spacecraft state is represented in Cartesian coordinates (Position and Velocity).

The algorithm proceeds as follows:
1.  **Jacobian Computation**: It computes the Jacobian matrix `J` representing the partial derivatives of the provided dispersion parameters with respect to the Cartesian state elements (x, y, z, vx, vy, vz).
2.  **Covariance Transformation**: It constructs a diagonal covariance matrix `P` in the parameter space (assuming input dispersions are independent in that space). It then transforms this into the Cartesian covariance matrix `C` using the linear mapping approximation:
`C = J_inv * P * J_inv^T`
where `J_inv` is the Moore-Penrose pseudo-inverse of `J`.
3.  **Decomposition**: It performs a Singular Value Decomposition (SVD) on `C` (or the full state covariance including mass, Cr, Cd) to obtain the square root of the covariance matrix, denoted as `L`.
`C = U * S * V^T = (V * sqrt(S)) * (V * sqrt(S))^T` implies `L = V * sqrt(S)`
4.  **Sampling**: To generate a sample state `X`, it draws a vector `Z` of independent standard normal variables (N(0, 1)) and applies the transformation:
`X = mu + L * Z`

# Correctness vs. Independent Sampling

One might ask: "Why not just sample each Cartesian coordinate independently using a normal distribution?"

1.  **Correlations**: Independent sampling of Cartesian coordinates assumes a diagonal covariance matrix, implying no correlation between position and velocity components. In orbital mechanics, states are highly correlated (e.g., velocity magnitude and radial distance are coupled by energy). `MvnSpacecraft` preserves these physical correlations by mapping the physically meaningful uncertainties (e.g., in SMA or Inclination) into the Cartesian space.
2.  **Geometry**: Uncertainties defined in orbital elements form complex shapes (like "bananas") in Cartesian space. A multivariate normal approximation in Cartesian space captures the principal axes and orientation of this uncertainty volume, which an axis-aligned bounding box (implied by independent sampling) effectively destroys.
3.  **Consistency**: By using the Jacobian transformation, we ensure that the generated samples, when mapped back to the parameter space (linearized), reproduce the input statistics (mean and standard deviation) provided by the user."""
    dispersions: typing.Any

    @staticmethod
    def from_spacecraft_cov(template: typing.Any, cov: typing.Any, mean: typing.Any) -> typing.Any:...

    def sample(self, count: typing.Any, seed: typing.Any=None) -> typing.Any:
        """Samples the multivariate distribution to generate a list of spacecraft states (up to 100k).

The Pseudo-Random Number Generator (PRNG) used is the Permuted Congruential Generator (PCG).
PCG is an excellent choice for Monte Carlo simulations because:
1. **Statistical Quality**: It passes difficult statistical tests (like TestU01), ensuring the generated numbers are random enough for high-fidelity simulations.
2. **Performance**: It is very fast and efficient, which is crucial when generating a large number of samples.
3. **Small State**: It has a small state size and is easy to seed, making it ideal for reproducible simulations.
4. **Reproducibility**: By providing a seed, the exact same sequence of spacecraft states can be generated, allowing for debugging and validation of Monte Carlo runs."""

    @staticmethod
    def zero_mean(template: typing.Any, dispersions: typing.Any) -> typing.Any:...

@typing.final
class StateDispersion:
    """A dispersions configuration, allows specifying min/max bounds (by default, they are not set)"""
    mean: typing.Any
    param: typing.Any
    std_dev: typing.Any

    @staticmethod
    def zero_mean(param: typing.Any, std_dev: typing.Any) -> typing.Any:...

@typing.final
class StateParameter:
    BLTOF: type = ...
    BdotR: type = ...
    BdotT: type = ...
    Cd: type = ...
    Cr: type = ...
    DryMass: type = ...
    Element: type = ...
    Epoch: type = ...
    GuidanceMode: type = ...
    Isp: type = ...
    PropMass: type = ...
    Thrust: type = ...
    ThrustInPlane: type = ...
    ThrustOutOfPlane: type = ...
    ThrustX: type = ...
    ThrustY: type = ...
    ThrustZ: type = ...
    TotalMass: type = ...