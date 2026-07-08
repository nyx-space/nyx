1. **Add `PositionDevice` class binding**: Edit `nyx-core/src/od/position/mod.rs` to add `pyclass(from_py_object)` macro.
2. **Add `PositionDevice` Python binding**: Create `nyx-core/src/od/position/python.rs` with `#[pymethods]` implementation. Add `pub mod python;` to `nyx-core/src/od/position/mod.rs`.
3. **Add `SpacecraftPositionKalmanOD`**: Add `pub type SpacecraftPositionKalmanOD = self::process::KalmanODProcess<SpacecraftDynamics, nalgebra::Const<1>, nalgebra::Const<3>, position::PositionDevice>;` to `nyx-core/src/od/mod.rs`.
4. **Export Position bindings**: Edit `nyx-py/src/lib.rs` to import and export `PositionDevice`, `PositionTrackingArcSim`, `PySpacecraftPositionODProcess`, and `PySpacecraftPositionODSolution`.
5. **Implement PySpacecraftPositionODProcess and PySpacecraftPositionODSolution**: Add these in `nyx-py/src/py_od.rs`, mirroring the existing GroundStation ones but for `PositionDevice`.
6. **Update `.pyi` type hints**: Add the new classes and methods to `nyx-py/nyx_space/orbit_determination.pyi`.
7. **Complete pre commit steps**: Ensure proper testing, verification, review, and reflection are done before submitting.
8. **Submit the change**.
