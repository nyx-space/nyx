use crate::celestia::{bodies, Cosm, Orbit};
use crate::dimensions::Vector3;
use crate::dynamics::orbital::OrbitalDynamics;
use crate::errors::NyxError;
use crate::propagators::{PropOpts, Propagator};
use crate::time::{Duration, TimeUnit};
use crate::tools::lambert::{standard, TransferKind};
use std::sync::mpsc::channel;

#[derive(Copy, Clone, Debug)]
pub struct Node {
    pub orbit: Orbit,
    pub ctrl: Vector3<f64>,
}

/// This trait determines the initial control nodes for direct multiple shooting.
pub trait Heuristic {
    fn nodes(&self, start: Orbit, end: Orbit) -> Result<Vec<Node>, NyxError>;
}

/// Lambert heuristic uses a Lambert transfer to find the shortest path between two nodes.
/// TODO: It will then use mesh refinement combined with a z-score to adequately place the control nodes.
/// NOTE: This is very likely a bad heuristic. The propagator tolerance is to the meter accuracy.
pub struct LambertHeuristic {
    pub tof: Duration,
}

impl Heuristic for LambertHeuristic {
    fn nodes(&self, start: Orbit, end: Orbit) -> Result<Vec<Node>, NyxError> {
        // Start by solving the Lambert solution between these nodes.
        let sol = standard(
            start.radius(),
            end.radius(),
            self.tof.in_seconds(),
            start.frame.gm(),
            TransferKind::Auto,
        )?;

        // Get into the transfer orbit.
        let mut start_tf = start;
        start_tf.vx = sol.v_init[0];
        start_tf.vy = sol.v_init[1];
        start_tf.vz = sol.v_init[2];

        // Propagate for the TOF
        let cosm = Cosm::de438();
        let bodies = vec![bodies::EARTH_MOON, bodies::SSB, bodies::JUPITER_BARYCENTER];
        let dynamics = OrbitalDynamics::point_masses(start.frame, &bodies, &cosm);

        // Create the channel to receive all of the details.
        let (tx, rx) = channel();
        let prop_setup = Propagator::rk89(&dynamics, PropOpts::with_tolerance(1e-3));
        let mut prop = prop_setup.with(start_tf);
        prop.tx_chan = Some(tx);
        prop.until_time_elapsed(self.tof)?;

        let mut node_states = Vec::new();
        while let Ok(state) = rx.try_recv() {
            node_states.push(state);
        }

        // Now that we have the nodes, let's compute the direction of each control
        let mut nodes = Vec::new();
        for i in 0..node_states.len() - 1 {
            let this = &node_states[i];
            let next = &node_states[i + 1];
            let delta = next.radius() - this.radius();
            let ctrl = delta / delta.norm();
            nodes.push(Node { orbit: *this, ctrl });
        }

        Ok(nodes)
    }
}

#[test]
fn lambert_heuristic() {
    use crate::time::Epoch;
    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);
    let dest_dt = Epoch::from_gregorian_utc_at_noon(2020, 1, 1);
    let delta_t = dest_dt - orig_dt;
    println!("{}", delta_t);

    let h = LambertHeuristic { tof: delta_t };

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let end = Orbit::keplerian(24_000.0, 0.1, 36.0, 60.0, 60.0, 180.0, dest_dt, eme2k);

    let start = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);

    let nodes = h.nodes(start, end).unwrap();

    println!("{}", nodes.len());
}
