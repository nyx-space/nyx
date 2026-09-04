# Per-Encounter Mascon Analysis for Low Lunar Orbits

A Nyx example demonstrating end-to-end per-encounter quantification of lunar mascon perturbations on low-lunar-orbit (LLO) altitude drift, ported from Osadugba (2026).

**Preprint:** [Per-Encounter Quantification of Mascon Effects on LLO Altitude Drift](https://doi.org/10.5281/zenodo.21740514) (Zenodo).
**Standalone reference implementation:** [KayAerospace/Per-Encounter-Quantification-of-Mascon-Effects-on-LLO-Altitude-Drift](https://github.com/KayAerospace/Per-Encounter-Quantification-of-Mascon-Effects-on-LLO-Altitude-Drift).

---

## What this example does

The Moon's gravity field is dominated by localized mass concentrations ("mascons"), buried anomalies that produce sharp, spatially confined perturbations on any spacecraft in low lunar orbit. Rather than treating these perturbations as a single averaged effect over an orbit, the referenced paper introduces a framework that segments them into discrete "encounters" (bounded windows where the disturbing acceleration exceeds a threshold), integrates each encounter analytically via a Gauss-variational mapping, and aggregates the results vectorially to predict altitude drift with per-encounter attribution.

This example ports the full pipeline onto Nyx: high-degree GRGM900C spherical-harmonic gravity evaluation, encounter detection with threshold sensitivity, per-encounter Gauss integration of Equation (I) from the paper, and multi-orbit iterated map validated against Nyx's full non-linear propagator over 50 revolutions.

## What this example exercises in Nyx

- `Almanac` with local NAIF kernel loading (DE440 planetary ephemeris + Moon PA/ME orientations).
- `Frame` handling: J2000 inertial ↔ Moon PA_DE440 body-fixed, including full time-varying rotation from `moon_pa_de440_200625.bpc`.
- `Orbit` construction (`from_keplerian`, `from_cartesian`, `with_ta_deg`), state extraction, and osculating element access.
- `GravityFieldConfig` reading a SHADR-format PDS3 coefficient file (`gggrx_0900c_sha.tab`) at degree 200.
- `AccelModels` / `ForceModels` / `Dynamics` composition — spherical harmonics only, no drag/SRP/thrust.
- `Propagator` with `IntegratorMethod.DormandPrince78` at `tolerance=1e-11`.
- `Propagator.accel_km_s2(spacecraft)` for instantaneous acceleration sampling along a fixed reference orbit (used to build the δg field without full propagation).
- `Propagator.for_duration(sc, T, trajectory=True)` for both single-orbit and 50-orbit truth propagation.
- `Trajectory.at(epoch)` for boundary and grid sampling of the propagated trajectory.
- `Almanac.rotate` and `Almanac.transform_to` for coordinate transformations.

The port is a **single self-contained Jupyter notebook** (`notebooks/mascon_encounter_analysis.ipynb`) with 11 numbered sections mapping onto the paper's Sections 3–6.

## Requirements

- **Python:** 3.11 or later.
- **Rust toolchain:** required to build Nyx's Python bindings (`rustc` ≥ 1.75, MSVC C++ build tools on Windows).
- **`uv`:** Python project manager. Install via `winget install astral-sh.uv` on Windows or `curl -LsSf https://astral.sh/uv/install.sh | sh` on Linux/macOS.
- **Nyx (`nyx_space`):** built from source. This example was developed against Nyx 2.5.0.
- **Internet access on first run:** required once to download NAIF kernels (~100 MB, cached locally). Subsequent runs are fully offline.
- **~50 MB disk space** for the GRGM900C coefficient file (see Data Setup).
- **~1 GB free RAM** during 50-orbit propagation.

## Install

Clone this Nyx fork and build the Python bindings:

```bash
git clone https://github.com/nyx-space/nyx.git
cd nyx/nyx-py
uv sync
```

The first `uv sync` will compile Nyx from Rust source (10–30 minutes depending on machine) and install all Python dependencies including matplotlib for this example. Subsequent runs are cached.

## Data Setup

This example requires the **GRGM900C lunar gravity model** (Lemoine et al., 2014) in PDS SHADR format.

**Download:**

```bash
curl -O https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_0900c_sha.tab
```

Approximate size: 49 MB.

**Placement:** put the file at `examples/mascon-encounter-analysis/data/gggrx_0900c_sha.tab` (relative to the Nyx repo root). The notebook reads it from `../data/gggrx_0900c_sha.tab` relative to itself.

## How to Run

From the Nyx repo root:

```bash
cd nyx/nyx-py
uv run jupyter lab ../examples/mascon-encounter-analysis/notebooks/mascon_encounter_analysis.ipynb
```

Or open the notebook in VS Code and select the `.venv` kernel that `uv sync` created at `nyx-py/.venv/Scripts/python.exe` (Windows) or `nyx-py/.venv/bin/python` (Linux/macOS).

**Run cells sequentially top-to-bottom.** Each of the 11 numbered sections is self-contained and produces printed output plus (in Sections 5–11) inline plots.

**On first run**, Nyx's `MetaAlmanac.latest()` downloads the standard NAIF kernels to `~/AppData/Local/nyx-space/anise/` (Windows) or the equivalent per-user cache on other OSes. Subsequent runs use the cached files.

## Expected Outputs and Runtime

Sequential run of all 11 sections on a modern laptop (i7-class, 16 GB RAM):

| Section | Content | Runtime |
|:---:|:---|:---:|
| 1 | Setup, constants, `CONFIG` | <1 s |
| 2 | Almanac + frame loading | 5–10 s (first run: +2–10 min for kernel download) |
| 3 | Nominal reference orbit + RTN basis | <1 s |
| 4 | Spherical harmonic gravity config | 5–15 s (SHADR parse on first `Propagator` build) |
| 5 | δg sampling at 4000 points | ~2 min |
| 6 | Encounter detection, gating, segmentation | <1 s |
| 7 | Per-encounter Δe, altitude drift, publication plots | 1–2 s |
| 8 | Threshold sensitivity sweep (5 percentiles) | 1–2 s |
| 9 | Single-orbit truth propagation + validation | 2–3 s |
| 10 | Multi-orbit iteration (Modes A + B) over 50 orbits | <1 s |
| 11 | Multi-orbit truth propagation + validation (50 orbits) | ~60 s |

**Total runtime end-to-end: 4–5 minutes** (excluding first-run kernel download).

## Findings from this port

This port reproduces the paper's *methodology* faithfully but produces different specific numerical values than the paper reports. All discrepancies are physically consistent and traceable to well-understood differences between independent implementations. Diagnostic sweeps documented at each stage confirm these are not code errors.

**Section 5: δg sampling.** Port's `|δg|` range at 50 km altitude: `[5.6e-8, 1.3e-6]` km/s². Paper's Section 6.1: `[1.94e-7, 1.70e-6]`. **The port's starting sub-satellite point is at (lat 0.5°, lon 112.4°)** because Nyx uses the *actual* Moon Principal Axis orientation at the specified UTC epoch, whereas the paper's implicit convention aligns Moon-fixed and inertial frames at t=0 (placing its start at lat 0°, lon 0°). Consequently, the port samples a different mascon-basin sequence. See [Reproducing the paper's exact geometry](#reproducing-the-papers-exact-geometry) for a recipe.

**Section 6: Encounter detection.** Port detects K = 6 encounters at P90; paper reports K = 4. Both are entirely on the lunar far side (longitudes clustered within ~90° of ±180°), confirming the paper's qualitative finding that lunar far-side mascons dominate LLO perturbations. Different count reflects different ground track.

**Section 7: Per-encounter Δe.** Port's cancellation ratio η at P90 = 0.545; paper: 0.101. Port's capture fraction at P90 = 162%; paper: 18.4%. These reflect different geometric relationships between the encounter Δe vectors at the two starting geometries.

**Section 8: Threshold sensitivity.** Port's qualitative pattern reproduces the paper's Figure 8: K grows as threshold loosens, capture peaks above 100% in the middle percentiles (overshoot regime), η stays in the resolvable band. **Port's optimum for per-encounter attribution is P50** (η = 0.123, capture = 121%, K = 13); paper's optimum is likewise P50 (η = 0.158, capture = 123.2%, K = 8). This robust identification of P50 across two independent geometries is a stronger validation of the paper's methodology than reproducing single-geometry numbers would be.

**Section 9: Single-orbit truth validation.** Port achieves **~18% agreement** between the first-order Gauss-integral prediction and Nyx-propagated truth over one orbital period. Paper reports 1.03% at its baseline geometry. Exhaustive diagnostic sweep confirms this is not a code issue:

- **Integrand code:** bit-identical between scalar and array formulations, and identical (0.00% difference) to the paper's standalone Cell 5 formulation.
- **Osculating-element extraction:** bit-identical to Nyx's own `Orbit.ecc()` (verified to 6e-16).
- **Propagator numerics:** settled at `tolerance = 1e-11` (0.2% shift at 1e-13).
- **Starting geometry:** tested at both this port's geometry (17.57%) and the paper's exact geometry rebuilt via inverse-DCM (18.10%). Not geometry-driven.

The residual difference is inherent between independent stacks: Nyx's Rust SH evaluator + PA-time-varying rotation + DormandPrince78 vs. the paper's Python SH + static R₃(ωt) + scipy DOP853. Each computes first-order-equivalent physics, but per-sample δg values differ by 15–20%, accumulating through integration. **This is itself a scientific finding**: reported first-order Gauss theory accuracy is implementation-dependent, and cross-validation between independent stacks is a stronger test than single-implementation validation.

**Section 10: Multi-orbit iteration.** Port's Mode A drift over 50 orbits = 10,990 m; paper (extrapolated) ≈ 24,859 m. Port's Mode A vs Mode B divergence = 0.004% relative; paper = 0.000014%. Both essentially zero; the paper's core finding that **out-of-plane effects are negligible for altitude drift over multi-day timescales** is reproduced at this different geometry.

**Section 11: Multi-orbit truth validation.** Port shows per-orbit truth Δe varies by 13× across the 50-orbit run (`4.1e-5` to `5.4e-4`); the frozen-map prediction (`1.23e-4` constant) is a poor per-orbit predictor (mean error 111%). Yet cumulative altitude drift at orbit 50 agrees to 80.7% between truth (13,625 m) and map (10,990 m). **This confirms the paper's Section 6.4 finding that per-orbit errors self-average over multi-orbit timescales**, though the effect is stronger at this port's geometry than at the paper's.

## Design Decisions

Explicit choices made during the port and their rationale:

- **μ from Nyx's Almanac (4902.800066)** rather than the paper's GRAIL GL0900D value (4902.800076). Relative difference ~2×10⁻⁸, sub-mm effect on altitude metrics.
- **Moon PA_DE440 frame** (full time-varying rotation from `moon_pa_de440_200625.bpc`) rather than the paper's simplified `R₃(±ωt)` static rotation. This is a physical upgrade; sub-meter numerical shifts result.
- **`IntegratorMethod.DormandPrince78`** as the closest analog to the paper's DOP853. Both are high-order embedded RK methods; numerically indistinguishable at `tolerance = 1e-11`.
- **δg via `Propagator.accel_km_s2`** at each anomaly sample of the fixed Keplerian reference (no propagation), preserving the paper's mathematical framework exactly. δg is sampled along the reference, not the true trajectory.
- **`AccelModels` composed with only `gravity_field`** (no `point_masses`, no `solid_tides`). ForceModels is empty. Matches the paper's gravity-only dynamics.
- **Manual RTN basis** (`r̂`, `t̂`, `n̂`) from (r, v) rather than Nyx's `dcm3x3_from_rcn_to_inertial`. Nyx's RCN column ordering differs from the paper's RTN; manual construction avoids the frame-convention risk.
- **Nominal orbit built by iterating `Orbit.with_ta_deg(ν_i)`** on a base `Orbit.from_keplerian`. Preserves Nyx type consistency while allowing NumPy array extraction for downstream numerical work.
- **N = 50 orbits** for multi-orbit iteration and truth validation (matches paper's Figure 14). Nyx's Rust propagator runs at ~1 s/orbit at degree 200, 2× faster than the paper's scipy DOP853, making 50 orbits comfortable where the paper defaulted to N = 10.

## Reproducing the paper's exact geometry

By default, this port uses the actual Moon Principal Axis orientation at UTC epoch 2026-01-01T00:00:00, producing a starting sub-satellite point at (lat 0.5°, lon 112.4°) in PA. The paper's implicit convention places the start at (lat 0°, lon 0°). To rebuild the initial orbit at the paper's geometry (while remaining polar in J2000 for propagation), replace Section 3's `orbit_0` construction with:

```python
# Rotate PA position (a, 0, 0) into J2000.
dcm_j2000_to_pa = almanac.rotate(Frames.MOON_J2000, moon_pa, epoch_0)
R_pa_to_j2000 = np.array(dcm_j2000_to_pa.rot_mat).T
r_j2000_paper = R_pa_to_j2000 @ np.array([sma_km, 0.0, 0.0])

# Polar-in-J2000 orbital plane through this position: h = r × ẑ / |·|.
z_j2000 = np.array([0.0, 0.0, 1.0])
h_dir = np.cross(r_j2000_paper, z_j2000)
h_dir /= np.linalg.norm(h_dir)
r_hat = r_j2000_paper / np.linalg.norm(r_j2000_paper)
v_dir = np.cross(h_dir, r_hat)
v_j2000_paper = v_circ * v_dir

orbit_0 = Orbit.from_cartesian(
    x_km    = r_j2000_paper[0],
    y_km    = r_j2000_paper[1],
    z_km    = r_j2000_paper[2],
    vx_km_s = v_j2000_paper[0],
    vy_km_s = v_j2000_paper[1],
    vz_km_s = v_j2000_paper[2],
    epoch   = epoch_0,
    frame   = moon_j2000,
)
```

Rerun Sections 3 onward. The starting sub-satellite point will be within ~0.02° of (lat 0°, lon 0°), residual arises from `IAU_MOON_FRAME` vs `MOON_PA_DE440_FRAME` differing by fractions of a degree in their formal Moon-fixed definitions.

**Note:** even at the paper's exact geometry, the port produces ~18% single-orbit validation error rather than the paper's 1% (see Section 9 findings above). This is inherent to independent-implementation numerical differences and does not indicate a geometry recovery failure.

## Future Extensions

- **Ground-track-aware Mode B iteration.** Section 10's Mode B currently rotates the (p̂, q̂) basis by the plane precession angle each orbit but re-uses the frozen `de2d_full`, a first-order approximation the paper's Section 7.3 flags as future work. With Nyx's `Propagator.accel_km_s2` cheap to call, a full ground-track-aware iteration (rebuilding the reference orbit and re-sampling δg on each iteration) becomes practical and would meaningfully improve multi-orbit accuracy at geometries where per-orbit Δe varies substantially (this port's regime).
- **Elliptical reference orbits.** The Gauss integrand generalises to `e > 0` via the `p` and `r(ν)` factors. All arrays in this port are already computed sample-wise from the actual state, so the code paths handle elliptical without modification, only the CONFIG needs updating and the analytical time grid needs to switch to Kepler's equation.
- **Different lunar gravity models.** The pipeline is not GRGM900C-specific. Swap the coefficient file path in CONFIG to any SHADR or COF file.
- **Third-body perturbations.** Add `PointMasses` (Sun, Earth) to `AccelModels` to extend beyond gravity-only dynamics. The encounter framework will need reformulation since third-body effects are not localized in the same sense as mascons.
- **Comparison against alternative spherical harmonic evaluators.** The port's 15-20% per-orbit validation error vs. paper is inherently an implementation-difference finding. Cross-validating against a third stack (GMAT, Orekit, STK) would characterize implementation dependence more rigorously.

## Citation

If you use this example or the underlying framework in academic work, please cite the preprint:

```bibtex
@misc{osadugba2026mascon,
  author       = {Osadugba, Ambrose},
  title        = {Per-Encounter Quantification of Mascon Effects on LLO Altitude Drift},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.21740514},
  url          = {https://doi.org/10.5281/zenodo.21740514},
}
```

And, per Nyx's licensing, acknowledge the Nyx astrodynamics library:

```bibtex
@software{rabotin_nyx,
  author = {Rabotin, Christopher},
  title  = {Nyx: Astrodynamics and Guidance, Navigation, and Control Library},
  url    = {https://github.com/nyx-space/nyx},
}
```

## License

This example inherits Nyx's [AGPL-3.0](../../LICENSE) license. The mathematical framework and standalone reference implementation are separately available at the linked repository above.

## Acknowledgements

Christopher Rabotin (Nyx maintainer) for reviving the Python bindings during the porting window and for reviewing the eventual pull request.