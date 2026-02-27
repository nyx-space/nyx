# Track Specification: Range Discrepancy Fix

## Problem
There is a reported range calculation discrepancy between `InterlinkTxSpacecraft` and `GroundStation` tracking when measuring an asset (e.g., a rover) on the surface of a celestial body. The user suspects that `InterlinkTxSpacecraft` provides an incorrect range when light-time (aberration) correction is enabled.

## Root Cause Hypothesis
`InterlinkTxSpacecraft::measure_instantaneous` computes the relative position vector as:
`let rho_tx_frame = rx_in_tx_frame.radius_km - observer.orbit.radius_km;` where `rx_in_tx_frame` is obtained via `almanac.transform_to(rx.orbit, observer.orbit.frame, self.ab_corr)`.

If `ab_corr` is `Aberration::LT`, `transform_to` corrects the light-time from the receiver to the *origin* of `observer.orbit.frame` (typically a celestial body center), not to the transmitter's position. This leads to an incorrect relative position vector when the transmitter is not at the origin of its own frame.

## Proposed Fix
The implementation should correctly compute the relative state of the receiver with respect to the transmitter, accounting for light-time between the two specific objects, rather than between the receiver and the frame origin.

## Verification Plan
1.  **Reproduction Test:** Create a test case mimicking the user's scenario: an orbiter and a rover on a celestial body. Use a body-fixed frame for the rover as described. Compare the range from an `InterlinkTxSpacecraft` (orbiter to rover) and a `GroundStation` (rover to orbiter).
2.  **Validation:** Ensure that the range and doppler values match between both implementations when configured identically.
