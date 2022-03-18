/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

extern crate serde;
extern crate serde_derive;

use self::serde_derive::{Deserialize, Serialize};
use crate::cosmic::{Cosm, Frame};
use crate::md::objective::Objective;
use crate::md::ui::StateParameter;
use crate::time::Epoch;
use crate::NyxError;
use std::convert::Into;
use std::str::FromStr;
use std::sync::Arc;

use super::multishoot::MultishootNode;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct NodesSerde {
    pub nodes: Vec<NodeSerde>,
}

impl NodesSerde {
    pub fn to_node_vec(&self, cosm: Arc<Cosm>) -> Result<Vec<Node>, NyxError> {
        let mut rtn = Vec::with_capacity(self.nodes.len());
        for n in &self.nodes {
            rtn.push(n.to_node(cosm.clone())?)
        }
        Ok(rtn)
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct NodeSerde {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vmag: Option<f64>,
    pub epoch: String,
    pub frame: String,
}

impl NodeSerde {
    pub fn to_node(&self, cosm: Arc<Cosm>) -> Result<Node, NyxError> {
        let frame = cosm.try_frame(self.frame.as_str())?;
        let epoch = Epoch::from_str(&self.epoch)?;

        Ok(Node {
            x: self.x,
            y: self.y,
            z: self.z,
            vmag: self.vmag.unwrap_or(0.0),
            frame,
            epoch,
        })
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Node {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vmag: f64,
    pub epoch: Epoch,
    pub frame: Frame,
}

impl Node {
    pub fn rmag(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }
}

impl MultishootNode<3> for Node {
    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn update_component(&mut self, component: usize, add_val: f64) {
        match component {
            0 => self.x += add_val,
            1 => self.y += add_val,
            2 => self.z += add_val,
            3 => self.vmag += add_val,
            _ => unreachable!(),
        }
    }
}

#[allow(clippy::from_over_into)]
impl Into<[Objective; 3]> for Node {
    fn into(self) -> [Objective; 3] {
        [
            Objective::new(StateParameter::X, self.x),
            Objective::new(StateParameter::Y, self.y),
            Objective::new(StateParameter::Z, self.z),
        ]
    }
}

impl MultishootNode<4> for Node {
    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn update_component(&mut self, component: usize, add_val: f64) {
        match component {
            0 => self.x += add_val,
            1 => self.y += add_val,
            2 => self.z += add_val,
            3 => self.vmag += add_val,
            _ => unreachable!(),
        }
    }
}

#[allow(clippy::from_over_into)]
impl Into<[Objective; 4]> for Node {
    fn into(self) -> [Objective; 4] {
        [
            Objective::new(StateParameter::X, self.x),
            Objective::new(StateParameter::Y, self.y),
            Objective::new(StateParameter::Z, self.z),
            Objective::new(StateParameter::Vmag, self.vmag),
        ]
    }
}

#[allow(clippy::from_over_into)]
impl Into<NodeSerde> for Node {
    fn into(self) -> NodeSerde {
        NodeSerde {
            x: self.x,
            y: self.y,
            z: self.z,
            vmag: Some(self.vmag),
            frame: self.frame.to_string(),
            epoch: self.epoch.to_string(),
        }
    }
}

#[allow(clippy::from_over_into)]
impl Into<NodeSerde> for &Node {
    fn into(self) -> NodeSerde {
        NodeSerde {
            x: self.x,
            y: self.y,
            z: self.z,
            vmag: Some(self.vmag),
            frame: self.frame.to_string(),
            epoch: self.epoch.to_string(),
        }
    }
}

#[test]
fn test_nodeserde() {
    extern crate toml;

    let str_nodes = r#"[[nodes]]
x = -394.37164017582654
y = -80.02184491079583
z = -1702.1160791417442
vmag = 0.0
epoch = "2023-11-25T14:11:46.789000034 UTC"
frame = "Moon J2000"

[[nodes]]
x = -381.68254116206856
y = -48.21573534985666
z = -1705.829637126235
vmag = 0.0
epoch = "2023-11-25T14:12:06.789000034 UTC"
frame = "Moon J2000"

[[nodes]]
x = -368.8474537620047
y = -16.401929604226403
z = -1708.8692139449731
vmag = 0.0
epoch = "2023-11-25T14:12:26.789000034 UTC"
frame = "Moon J2000"
"#;

    let toml_nodes: NodesSerde = toml::from_str(str_nodes).unwrap();

    let cosm = Cosm::de438();

    let nodes = toml_nodes.to_node_vec(cosm).unwrap();

    dbg!(&nodes);

    let v = NodesSerde {
        nodes: nodes.iter().map(|n| n.into()).collect::<Vec<NodeSerde>>(),
    };
    let toml_ser = toml::to_string(&v).unwrap();

    println!("GOT\n{}\n\nWANTED:{}", toml_ser, str_nodes);

    assert_eq!(toml_ser, str_nodes);
}
