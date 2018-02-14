extern crate nalgebra as na;
extern crate nyx;
use self::na::DVector;

fn two_body_dynamics(_t: f64, state: DVector<f64>) -> DVector<f64> {
    let radius = state.slice((0, 0), (3, 1)); // TODO: Change to compile time slice
    let velocity = state.slice((3, 0), (3, 1));
    let body_acceleration = (-398600.4 / radius.norm().powi(3)) * radius;
    DVector::from_iterator(6, velocity.iter().chain(body_acceleration.iter()).cloned())
}

#[test]
fn geo_day_prop() {
    extern crate nyx;
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::DVector;
    //let init_state = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);
    let init_state =
        DVector::from_row_slice(6, &[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
    let method = RKF54::from_options(Options::with_fixed_step(1.0));
    let mut prop = Propagator::new::<RKF54>(0.0, init_state.clone(), method);
    let state = prop.prop(two_body_dynamics);
    println!("delta state: \n {}", state - init_state);
    /*while prop.latest_time() < 3600.0 * 24.0 {
        let state = prop.prop(two_body_dynamics);
        println!("final state: \n {}", state);
    }*/

    /*
virtObj := CelestialObject{"virtObj", 6378.145, 149598023, 398600.4, 23.4, 0.00005, 924645.0, 0.00108248, -2.5324e-6, -1.6204e-6, 0, nil}
orbit := NewOrbitFromRV([]float64{-2436.45, -2436.45, 6891.037}, []float64{5.088611, -5.088611, 0}, virtObj)
startDT := time.Date(2017, 1, 1, 0, 0, 0, 0, time.UTC)
endDT := startDT.Add(24 * time.Hour).Add(time.Second)
NewPreciseMission(NewEmptySC("est", 0), orbit, startDT, endDT, Perturbations{}, time.Second, false, ExportConfig{}).Propagate()
expR := []float64{-5971.19544867343, 3945.58315019255, 2864.53021742433}
expV := []float64{0.049002818030, -4.185030861883, 5.848985672439}
if !floats.EqualApprox(orbit.rVec, expR, 1e-8) {
    t.Fatalf("Incorrect R:\ngot: %+v\nexp: %+v", orbit.rVec, expR)
}
if !floats.EqualApprox(orbit.vVec, expV, 1e-8) {
    t.Fatalf("Incorrect R:\ngot: %+v\nexp: %+v", orbit.vVec, expV)
}

R := []float64{f[0], f[1], f[2]}
V := []float64{f[3], f[4], f[5]}
tmpOrbit = NewOrbitFromRV(R, V, a.Orbit.Origin)
bodyAcc := -tmpOrbit.Origin.μ / math.Pow(Norm(R), 3)
_, _, i, Ω, _, _, _, _, u := tmpOrbit.Elements()
// Check if any impulse burn, and execute them if needed.
if maneuver, exists := a.Vehicle.Maneuvers[a.CurrentDT.Truncate(a.step)]; exists {
    if !maneuver.done {
        a.Vehicle.logger.Log("level", "info", "subsys", "astro", "date", a.CurrentDT, "thrust", "impulse", "v(km/s)", maneuver.Δv())
        Δv[0] += maneuver.R
        Δv[1] += maneuver.N
        Δv[2] += maneuver.C
        maneuver.done = true
        a.Vehicle.Maneuvers[a.CurrentDT.Truncate(a.step)] = maneuver
    }
}
Δv = Rot313Vec(-u, -i, -Ω, Δv)
// d\vec{R}/dt
fDot[0] = f[3]
fDot[1] = f[4]
fDot[2] = f[5]
// d\vec{V}/dt
fDot[3] = bodyAcc*f[0] + Δv[0]
fDot[4] = bodyAcc*f[1] + Δv[1]
fDot[5] = bodyAcc*f[2] + Δv[2]
*/
}
