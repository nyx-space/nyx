use arrow::array::{Array, LargeStringArray, StringArray};
use std::sync::Arc;

pub enum AbstractStringArray<'a> {
    Small(&'a StringArray),
    Large(&'a LargeStringArray),
}

impl<'a> AbstractStringArray<'a> {
    pub fn try_from(array: &'a Arc<dyn Array>) -> Option<Self> {
        if let Some(s) = array.as_any().downcast_ref::<StringArray>() {
            Some(AbstractStringArray::Small(s))
        } else {
            array
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .map(AbstractStringArray::Large)
        }
    }

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
