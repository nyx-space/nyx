/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::io::ConfigError;
use pyo3::types::{PyAny, PyDict, PyList};
use serde_yaml::{Mapping, Value};

/// Try to convert the provided PyAny into a SerDe YAML Value
pub fn pyany_to_value(any: &PyAny) -> Result<Value, ConfigError> {
    if let Ok(as_bool) = any.extract::<bool>() {
        Ok(Value::Bool(as_bool))
    } else if let Ok(as_str) = any.extract::<String>() {
        Ok(Value::String(as_str))
    } else if let Ok(as_f64) = any.extract::<f64>() {
        Ok(Value::Number(as_f64.into()))
    } else if let Ok(as_int) = any.extract::<i64>() {
        Ok(Value::Number(as_int.into()))
    } else if let Ok(as_int) = any.extract::<i32>() {
        Ok(Value::Number(as_int.into()))
    } else if let Ok(as_list) = any.downcast::<PyList>() {
        let mut seq = Vec::new();
        for item in as_list.iter() {
            seq.push(pyany_to_value(item)?);
        }
        Ok(Value::Sequence(seq))
    } else if let Ok(as_dict) = any.downcast::<PyDict>() {
        let mut map = Mapping::new();
        for item_as_list in as_dict.items() {
            let k = pyany_to_value(item_as_list.get_item(0).or_else(|_| {
                Err(ConfigError::InvalidConfig {
                    msg: "could not get key from provided dictionary item".to_string(),
                })
            })?)?;
            let v = pyany_to_value(item_as_list.get_item(1).or_else(|_| {
                Err(ConfigError::InvalidConfig {
                    msg: "could not get value from provided dictionary item".to_string(),
                })
            })?)?;
            map.insert(k, v);
        }

        Ok(Value::Mapping(map))
    } else {
        Err(ConfigError::InvalidConfig {
            msg: "cannot convert input (not one of bool, int, float, str, List, Dict)".to_string(),
        })
    }
}
