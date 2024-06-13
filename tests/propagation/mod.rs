pub(crate) const GMAT_EARTH_GM: f64 = 398_600.441_5;
pub(crate) const GMAT_SUN_GM: f64 = 132_712_440_017.99;
pub(crate) const GMAT_MOON_GM: f64 = 4_902.800_582_147_8;

mod events;
mod propagators;
mod stm;
mod stopcond;
mod trajectory;
