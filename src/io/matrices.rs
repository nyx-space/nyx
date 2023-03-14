/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use either::Either;
use nalgebra::{Matrix2, Matrix6, Vector2, Vector6};
use serde::{Deserialize, Serialize};

// TODO: Use a macro to generate all of these implementations

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Matrix2Serde {
    #[serde(with = "either::serde_untagged")]
    inner: Either<Diag2, Mat2>,
}

impl Matrix2Serde {
    pub fn to_matrix(&self) -> Matrix2<f64> {
        match self.inner {
            Either::Left(diag) => Matrix2::from_diagonal(&Vector2::from_iterator(diag.0)),
            Either::Right(mat2) => {
                let mut flat: [f64; 4] = [0.0; 4];
                for (i, row) in mat2.0.iter().enumerate() {
                    for (j, val) in row.iter().enumerate() {
                        flat[4 * i + j] = *val;
                    }
                }
                Matrix2::from_row_slice(&flat)
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Diag2([f64; 2]);

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Mat2([[f64; 2]; 2]);

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Matrix6Serde {
    #[serde(with = "either::serde_untagged")]
    inner: Either<Diag6, Mat6>,
}

impl Matrix6Serde {
    pub fn to_matrix(&self) -> Matrix6<f64> {
        match self.inner {
            Either::Left(diag) => Matrix6::from_diagonal(&Vector6::from_iterator(diag.0)),
            Either::Right(mat6) => {
                let mut flat: [f64; 36] = [0.0; 36];
                for (i, row) in mat6.0.iter().enumerate() {
                    for (j, val) in row.iter().enumerate() {
                        flat[6 * i + j] = *val;
                    }
                }
                Matrix6::from_row_slice(&flat)
            }
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Diag6([f64; 6]);

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Mat6([[f64; 6]; 6]);

#[test]
fn test_serde2() {
    use serde_yaml;

    let m_diag = Matrix2Serde {
        inner: Either::Left(Diag2([1.0, 2.0])),
    };

    println!("Diag -- \n{}", serde_yaml::to_string(&m_diag).unwrap());
    // Load from one line list
    let diag_s = "[1.0, 2.0]";
    let diag_loaded: Matrix2Serde = serde_yaml::from_str(diag_s).unwrap();
    assert_eq!(diag_loaded, m_diag);

    let m_full = Matrix2Serde {
        inner: Either::Right(Mat2([[1.0, 2.0]; 2])),
    };

    // Serialization will print this as an exhaustive list of lists.
    println!("Full -- \n{}", serde_yaml::to_string(&m_full).unwrap());
    // Load from list
    let full_mat = r#"
- [1.0, 2.0] # Row 1
- [1.0, 2.0] # Row 2
    "#;

    let full_loaded: Matrix2Serde = serde_yaml::from_str(full_mat).unwrap();

    assert_eq!(full_loaded, m_full);
}

#[test]
fn test_serde6() {
    use serde_yaml;

    let m_diag = Matrix6Serde {
        inner: Either::Left(Diag6([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])),
    };

    println!("Diag -- \n{}", serde_yaml::to_string(&m_diag).unwrap());
    // Load from one line list
    let diag_s = "[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]";
    let diag_loaded: Matrix6Serde = serde_yaml::from_str(diag_s).unwrap();
    assert_eq!(diag_loaded, m_diag);

    let m_full = Matrix6Serde {
        inner: Either::Right(Mat6([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]; 6])),
    };

    // Serialization will print this as an exhaustive list of lists.
    println!("Full -- \n{}", serde_yaml::to_string(&m_full).unwrap());
    // Load from list
    let full_mat = r#"
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 1
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 2
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 3
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 4
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 5
- [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] # Row 6
    "#;

    let full_loaded: Matrix6Serde = serde_yaml::from_str(full_mat).unwrap();

    assert_eq!(full_loaded, m_full);
}
