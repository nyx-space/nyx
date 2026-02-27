# Implementation Plan - Range Discrepancy Fix

## Phase 1: Research and Reproduction
- [ ] Task: Create a reproduction script or test case that demonstrates the range discrepancy between `InterlinkTxSpacecraft` and `GroundStation`.
    - [ ] Setup a scenario with an orbiter and a rover on the Moon.
    - [ ] Configure `InterlinkTxSpacecraft` with `Aberration::LT`.
    - [ ] Configure `GroundStation` with `light_time_correction: true`.
    - [ ] Assert that the ranges differ significantly.

## Phase 2: Implementation of Fix
- [ ] Task: Modify `InterlinkTxSpacecraft::measure_instantaneous` to use a robust relative state calculation.
    - [ ] Investigate if `anise` provides a direct method for relative state with aberration correction between two orbits.
    - [ ] If not, implement the light-time correction manually or by using a temporary frame centered on the transmitter.
- [ ] Task: Ensure Doppler (range-rate) is also correctly updated for the new relative state.

## Phase 3: Verification and Cleanup
- [ ] Task: Verify the fix using the reproduction test case.
- [ ] Task: Check for regressions in other tracking examples.
- [ ] Task: Conductor - User Manual Verification 'Phase 3' (Protocol in workflow.md)
