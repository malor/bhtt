#[cfg(test)]
#[macro_use]
extern crate approx;

mod bin;
mod histogram;

pub use bin::Bin;
pub use histogram::Histogram;
