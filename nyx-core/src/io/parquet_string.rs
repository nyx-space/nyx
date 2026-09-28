use arrow::{
    array::{Array, LargeStringArray, StringArray},
    error::ArrowError,
};
use std::sync::Arc;

pub enum AbstractStringArray<'a> {
    Small(&'a StringArray),
    Large(&'a LargeStringArray),
}

impl<'a> AbstractStringArray<'a> {
    pub fn value(&self, i: usize) -> &str {
        match self {
            AbstractStringArray::Small(s) => s.value(i),
            AbstractStringArray::Large(s) => s.value(i),
        }
    }

    pub fn is_null(&self, i: usize) -> bool {
        match self {
            AbstractStringArray::Small(s) => s.is_null(i),
            AbstractStringArray::Large(s) => s.is_null(i),
        }
    }

    pub fn is_valid(&self, i: usize) -> bool {
        match self {
            AbstractStringArray::Small(s) => s.is_valid(i),
            AbstractStringArray::Large(s) => s.is_valid(i),
        }
    }
}

impl<'a> TryFrom<&'a Arc<dyn Array>> for AbstractStringArray<'a> {
    type Error = ArrowError;

    fn try_from(array: &'a Arc<dyn Array>) -> Result<Self, Self::Error> {
        if let Some(s) = array.as_any().downcast_ref::<StringArray>() {
            Ok(AbstractStringArray::Small(s))
        } else {
            match array.as_any().downcast_ref::<LargeStringArray>() {
                Some(downcasted) => Ok(Self::Large(downcasted)),
                None => Err(ArrowError::CastError(
                    "column is neither StringArray nor LargeStringArray".to_string(),
                )),
            }
        }
    }
}
