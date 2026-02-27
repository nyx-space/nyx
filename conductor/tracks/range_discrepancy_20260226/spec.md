# Track Specification: Range Discrepancy Fix

## Problem
There is a reported range and doppler calculation discrepancy between `InterlinkTxSpacecraft` and `GroundStation` tracking. This happens even when aberration (light-time) correction is disabled.

## Root Cause
`InterlinkTxSpacecraft::measure_instantaneous` computes the relative state using naive vector subtraction and an incorrect assumption:
```rust
let rho_tx_frame = rx_in_tx_frame.radius_km - observer.orbit.radius_km;

// Compute the range-rate \dot ρ. Note that rx_in_tx_frame is already the relative velocity of rx wrt tx!
let range_rate_km_s = rho_tx_frame.dot(&rx_in_tx_frame.velocity_km_s) / rho_tx_frame.norm();
```

The comment "Note that rx_in_tx_frame is already the relative velocity of rx wrt tx!" is **incorrect** unless the observer is stationary at the origin of `observer.orbit.frame`. For an orbiter in EME2000 or a body-fixed frame, `rx_in_tx_frame.velocity_km_s` is the receiver's velocity in that frame, not the relative velocity. The correct relative velocity is `rx_in_tx_frame.velocity_km_s - observer.orbit.velocity_km_s`.

Furthermore, relying on manual vector subtraction is less robust than using high-level geometric functions provided by `anise`.

## Proposed Fix
Modify `InterlinkTxSpacecraft::measure_instantaneous` to:
1.  Correctly compute the relative position and velocity vectors by subtracting the observer's state from the receiver's state (both transformed to the same frame).
2.  Investigate using `anise`'s `relative_state` or similar high-level functions to ensure consistency with `GroundStation` and other tracking devices.
3.  Ensure the fix works correctly both with and without aberration corrections.

## Verification Plan
1.  **Reproduction Test:** Create a test case with an orbiter and a "rover" (mimicked by a body-fixed trajectory). Compare range and doppler from `InterlinkTxSpacecraft` and `GroundStation`.
2.  **Validation:** Verify that range and doppler values match between both implementations within acceptable numerical precision.
