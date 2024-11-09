use super::Base;

pub struct Kilograms(Base);

super::unit!(Kilograms, "kgs");

pub struct Pounds(Base);

super::unit!(Pounds, "lbs");

super::linear_conversion!(2.20462_f64, Kilograms, Pounds);
super::linear_conversion!(0.453592_f64, Pounds, Kilograms);
