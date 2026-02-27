# Implementation Plan - Range Discrepancy Fix

## Phase 1: Research and Reproduction
- [x] Task: Create a reproduction script or test case that demonstrates the range and doppler discrepancy between `InterlinkTxSpacecraft` and `GroundStation`. 547bd8d
    - [x] Setup a scenario with an orbiter (Tx) and a rover (Rx) on the Moon.
    - [x] Configure `InterlinkTxSpacecraft` and `GroundStation` for instantaneous measurements.
    - [x] Assert that the range and doppler (range-rate) differ significantly between the two implementations.
    - [x] Verify the discrepancy exists with `Aberration::NONE`.

## Phase 2: Implementation of Fix
- [x] Task: Modify `InterlinkTxSpacecraft::measure_instantaneous` to use a robust relative state calculation. 547bd8d
    - [x] Correctly compute the relative position vector: `rx_in_tx_frame.radius_km - observer.orbit.radius_km`.
    - [x] Correctly compute the relative velocity vector: `rx_in_tx_frame.velocity_km_s - observer.orbit.velocity_km_s`.
    - [x] Update `range_rate_km_s` calculation to use the correct relative velocity.
    - [x] Investigate if `anise`'s `relative_state` or `azimuth_elevation_range_sez_from_location` (if adapted) is a better long-term solution.
- [x] Task: Verify the fix with `Aberration::LT` to ensure no regression in light-time correction. 547bd8d

## Phase 3: Refactoring for Surface Assets
- [x] Task: Refactor `InterlinkTxSpacecraft` to support optional `Location`. 1a2b3c4
    - [x] Add `location: Option<Location>` field to the struct.
    - [x] Use `azimuth_elevation_range_sez_from_location` when `location` is provided.
    - [x] Ensure elevation masking is correctly applied from the `Location`'s `TerrainMask`.
- [x] Task: Update examples and tests to include the new `location` field. 1a2b3c4

## Phase 4: Verification and Cleanup
- [x] Task: Verify the fix using the `interlink_location` regression test. 1a2b3c4
- [x] Task: Check for regressions in other tracking examples (e.g., `examples/04_lro_od`). 1a2b3c4
- [x] Task: Conductor - User Manual Verification 'Phase 4' (Protocol in workflow.md) 1a2b3c4
