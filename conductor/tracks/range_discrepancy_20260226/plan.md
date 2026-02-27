# Implementation Plan - Range Discrepancy Fix

## Phase 1: Research and Reproduction
- [ ] Task: Create a reproduction script or test case that demonstrates the range and doppler discrepancy between `InterlinkTxSpacecraft` and `GroundStation`.
    - [ ] Setup a scenario with an orbiter (Tx) and a rover (Rx) on the Moon.
    - [ ] Configure `InterlinkTxSpacecraft` and `GroundStation` for instantaneous measurements.
    - [ ] Assert that the range and doppler (range-rate) differ significantly between the two implementations.
    - [ ] Verify the discrepancy exists with `Aberration::NONE`.

## Phase 2: Implementation of Fix
- [ ] Task: Modify `InterlinkTxSpacecraft::measure_instantaneous` to use a robust relative state calculation.
    - [ ] Correctly compute the relative position vector: `rx_in_tx_frame.radius_km - observer.orbit.radius_km`.
    - [ ] Correctly compute the relative velocity vector: `rx_in_tx_frame.velocity_km_s - observer.orbit.velocity_km_s`.
    - [ ] Update `range_rate_km_s` calculation to use the correct relative velocity.
    - [ ] Investigate if `anise`'s `relative_state` or `azimuth_elevation_range_sez_from_location` (if adapted) is a better long-term solution.
- [ ] Task: Verify the fix with `Aberration::LT` to ensure no regression in light-time correction.

## Phase 3: Verification and Cleanup
- [ ] Task: Verify the fix using the reproduction test case.
- [ ] Task: Check for regressions in other tracking examples (e.g., `examples/04_lro_od`).
- [ ] Task: Conductor - User Manual Verification 'Phase 3' (Protocol in workflow.md)
