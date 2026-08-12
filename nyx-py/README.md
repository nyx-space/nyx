# Nyx: High-Fidelity Flight Dynamics and Orbit Determination

_Blazing fast astrodynamics from mission concept to operations, analysis, and automation_

```sh
    pip install nyx_space
```

Requires Python 3.11 or later; [works on Linux, Windows, MacOS (Intel and arm64)](https://pypi.org/project/nyx_space/#files).

## High fidelity example

_Note:_ You must swap the paths in this example with data downloaded from the Github repo: <https://github.com/nyx-space/nyx/tree/master/data>.

```python
import logging

from nyx_space import DragData, ExportCfg, Mass, Spacecraft, SRPData
from nyx_space.anise import Aberration, MetaAlmanac
from nyx_space.anise.analysis import Event
from nyx_space.anise.astro import DataType, Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AccelModels,
    AtmDensity,
    Drag,
    Dynamics,
    ForceModels,
    GravityFieldConfig,
    IntegratorMethod,
    IntegratorOptions,
    Nrlmsise00Flags,
    PointMasses,
    Propagator,
    SolarPressure,
    SolidTides,
    SpaceWeatherData,
    StaticSpaceWeather,
)
from nyx_space.time import Duration, Epoch


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Load the universe via ANISE's MetaAlmanac (similar to SPICE's meta kernel)
    almanac = (
        MetaAlmanac("../data/02_config/ci_almanac.dhall")
        .process()
        .load("../data/01_planetary/earth_2025_250826_2125_predict.bpc")
    )

    # Acceleration models are independent of the vehicle mass.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[
                CelestialObjects.EARTH,
                CelestialObjects.MOON,
                CelestialObjects.SUN,
            ],
        ),
        gravity_field=GravityFieldConfig(
            degree=50,
            order=50,
            filepath="../data/01_planetary/EGM2008_to2190_TideFree_sha.gz",
            frame=Frames.EARTH_ITRF93.to_frameuid(),
        ),
        # IMPORTANT: Solid tides are defined in the body fixed frame, so we must initialize the model with these frames.
        solid_tides=SolidTides.earth_moon_system(
            Frames.IAU_EARTH_FRAME, moon_frame=Frames.IAU_MOON_FRAME, almanac=almanac
        ),
    )

    # Define Non-Gravitational Forces: force models depends on the vehicle mass.
    # Configure Solar Radiation Pressure targeting Earth and Moon frames, accounting for light-time aberration
    # because the direction of the photons depend on the position of the sun at time of emission.
    srp = SolarPressure(
        [Frames.EARTH_J2000, Frames.MOON_J2000], almanac, correction=Aberration("LT")
    )
    # Configure atmospheric drag using a standard exponential atmospheric model.
    # Load the Space Weather as provided by CelesTrak: https://celestrak.org/SpaceData/.
    # You should also provide a fallback for propagation that extends further from the
    # data set. This can be specified either as Solar Minimum, Solar Maximum, Solar Average,
    # or a custom number for the F10.7 (in SFU), Ap, and Kp values.
    # Data may be provided either as CSV or in non-archived gunzip (gz) format (decoded on the fly).
    weather = SpaceWeatherData(
        "../data/01_planetary/SpaceWeather-2021-01-01_2026-09-06.csv.gz",
        StaticSpaceWeather.SolarAverage(),
    )
    # NOTE While you must specify the flags for the NRLMSISE00 models, the default ones are typically what you want.
    # Call `help(Nrlmsise00Flag)` for details on available options.
    drag = Drag(
        AtmDensity.NRLMSISE00(weather, Nrlmsise00Flags()),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
    )

    # Combine accelerations and forces into a unified Dynamics model.
    dynamics = Dynamics(accel_models, force_models=ForceModels(srp, drag))

    # Initialize the Propagator
    # The Propagator encapsulates the dynamics, the math (RungeKutta89), and the universe (almanac).
    # It remains stateless regarding the spacecraft, allowing it to be reused.
    propagator = Propagator(
        dynamics, almanac, IntegratorMethod.RungeKutta89, IntegratorOptions()
    )

    # Define the Spacecraft Initial State
    # We define the starting Keplerian orbital elements.
    eme2k = almanac.frame_info(Frames.EME2000)
    orbit = Orbit.from_keplerian(
        sma_km=6800.0,
        ecc=0.0123,
        inc_deg=45.0,
        raan_deg=60.0,
        aop_deg=75.0,
        ta_deg=90.0,
        epoch=Epoch("2030-12-29 01:02:03 TDB"),
        frame=eme2k,
    )

    # Define the physical characteristics of the spacecraft necessary to compute drag and SRP.
    # Surfaces are in square meters, masses are in kilograms.
    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=123.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.0),
    )

    # Propagate until a specific orbital event (periapsis).
    # This demonstrates the root-finding capabilities of Nyx to stop exactly at an event,
    # bounded by a maximum duration to prevent infinite searching.
    result = propagator.until_event(
        my_sc, Event.periapsis(), max_duration=Duration("1 day")
    )
    logging.info(
        f"Periapsis found {result.state.orbit.epoch - orbit.epoch} after initial state."
    )

    # Export to SPICE BSP format
    ephem = result.trajectory.to_ephemeris("prop-to-peri", ExportCfg())
    ephem.write_spice_bsp(
        -200_000, "prop_to_peri.bsp", DataType.Type13HermiteUnequalStep
    )
```

More examples:

- <https://github.com/nyx-space/nyx/tree/master/nyx-py/examples>
- <https://github.com/nyx-space/nyx/tree/master/nyx-py/tests>
- [https://nyxspace.com/nyxspace/showcase/](https://nyxspace.com/nyxspace/showcase/?utm_source=pyreadme)

Flight dynamics operations demand an unforgiving synthesis of computational speed and physical fidelity.

While heritage tools provide the necessary precision for deep space and near-Earth navigation, their monolithic architectures and sluggish sequential processing bottleneck modern Monte Carlo analysis and autonomous operations.
Nyx answers this limitation by marrying a high-performance Rust computational engine with a productive Python interface, achieving sub-meter accuracy in high-fidelity scenarios. It is validated and proven to [match GMAT and industry-standard commercial astrodynamics suites to sub-meter tolerances](https://nyxspace.com/nyxspace/MathSpec/?utm_source=pyreadme).

## The Engine: Determinism and Precision

At its core, Nyx relies on rigorously validated physical models. Propagators support variable-step and fixed-step integrators (such as Runge-Kutta 8(9) and Dormand-Prince 7(8)), handling complex dynamics without wasting CPU cycles.

* **Gravitational Modeling:** Supports high-degree and order spherical harmonics (e.g., EGM2008 or any [SHADR file from NASA Planetary Data Service](https://pds-geosciences.wustl.edu/dataserv/gravity_models.htm)) and fully configurable solid tides compatible from the Cislunar systems and to gas giants via configurable Love numbers.
* **Non-Gravitational Forces:** Features robust Solar Radiation Pressure (SRP) and atmospheric drag models. Notably, the NRLMSISE-00 atmospheric model is fully compatible with the [CelesTrak (Cssi) SpaceWeather data](https://celestrak.org/SpaceData/).
* **Time and Frames:** Built on top of [`hifitime`](https://nyxspace.com/hifitime/?utm_source=pyreadme) and [`ANISE`](https://nyxspace.com/anise/?utm_source=pyreadme), Nyx guarantees strict nanosecond precision across 32,768 centuries, avoiding the floating-point inaccuracies common in legacy systems. It fluently handles GPST, ET/TDB, TAI, UTC, and [many more](https://docs.rs/hifitime/latest/hifitime/enum.TimeScale.html#variants).

## Orbit Determination: Statistical Rigor

The orbit determination module in Nyx allows estimation of position, velocity, coefficient of reflectivity, coefficient of drag, and ground station location and crust movement. Further enhancements are on the way.
The OD capabilities revolve around an Extended Kalman Filter (EKF), supporting both Deviation Tracking and Reference Update variants, coupled with backward smoothing.

* **Measurement Simulation:** Generate synthetic radiometric tracking data (Range, Doppler, Angles) across complex ground station networks. The simulator applies stochastic noise models, including white noise and Gauss-Markov bias processes, properly deriving thermal noise from Carrier-to-Noise density (C/N0) rather than artificially inflating velocity noise.
* **Consistency Validation:** A filter is useless if its confidence bounds do not reflect reality. Nyx natively implements Normalized Innovation Squared (NIS) and Normalized Estimation Error Squared (NEES) statistical tests, alongside Kolmogorov–Smirnov normality checks, ensuring covariance realism.
* **Interoperability:** Seamlessly ingest and export tracking data and trajectories using CCSDS standards (TDM and OEM) or Apache Parquet for high-throughput data pipelines and plotting.

## Monte Carlo and Concurrency: Bypassing the GIL

The true bottleneck of Python in flight dynamics is the Global Interpreter Lock (GIL). Nyx circumvents this by delegating heavy computation directly to Rust.

When executing a Monte Carlo dispersion, Nyx utilizes the `MvnSpacecraft` generator. Instead of naively sampling independent Cartesian coordinates, it computes the Jacobian matrix of the orbital elements, constructs a diagonal covariance matrix, and transforms this into Cartesian space via a Moore-Penrose pseudo-inverse and Singular Value Decomposition (SVD). This preserves the physical correlations of the state, correctly translating Keplerian uncertainty volumes into operational Cartesian bounds.

By calling `propagator.many_until_event(...)`, these dispersed states are propagated in parallel across all available CPU cores, completely bypassing the Python GIL for massive computational throughput.

In a benchmark stress-testing parallel execution, Nyx propagated an ensemble of 5,000 dispersed states for the James Webb Space Telescope across a 6.5-day high-fidelity deep-space trajectory, including gravitational pull from the Earth, Moon, Sun, and Jupiter barycenter, plus solar radiation pressure with Earth and Moon shadow modeling.
Nyx completed the full Monte Carlo run in 3 minutes 17 seconds on a 6-core Intel i7-8700K (12 threads) @ 3.70 GHz: this throughput is truly built for rapid, operational stochastic analysis rather than overnight batch jobs.

## Missions

Nyx has been used in various capacities on the following missions:

- [Firefly Blue Ghost series](https://fireflyspace.com/blue-ghost/), the only successful NASA CLPS mission
- [Masten XL-1](https://masten.aero/xelene-lunar-lander/), a former NASA CLPS mission
- [NASA CAPSTONE](https://www.nasa.gov/smallspacecraft/capstone/), a pioneer in spacecraft position, navigation, and timing

Nyx offers the operational robustness required to land a spacecraft on the Moon, providing FDOs with a toolchain that respects both the laws of physics and the value of engineering time.

## License

Nyx Space is provided under the AGPL-3.0 copyleft license.

Under the AGPL-3.0 terms, any user is permitted free and unlimited internal use of the library both for evaluation and operational deployment within an organization. In other words, there is no limitation to use Nyx as the core flight dynamics engine from preliminary design all the way to operations, provided your company maintains full internal control of the software stack and its users. The license's source-disclosure obligations only trigger if Nyx is integrated into a product or network service made accessible (directly or indirectly) to external third parties, e.g. if you build a commercial software package incorporating Nyx and distribute it externally (again, internal distribution does not trigger the AGPL propagation clause).
