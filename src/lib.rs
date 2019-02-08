#[macro_use]
extern crate failure;
extern crate libc;
extern crate sfcgal_sys;
#[macro_use]
extern crate enum_primitive_derive;
extern crate num_traits;
extern crate geo_types;
#[macro_use]
extern crate approx;

use sfcgal_sys::sfcgal_version;

mod conversion;
mod errors;
mod utils;
mod geometry;
mod coords;
pub use errors::Result;
pub use geometry::{SFCGeometry, GeomType};
pub use conversion::TryInto;
pub use coords::{CoordSeq};


/// Convert object to a SFCGAL geometry.
pub trait ToSfcgal {
    fn to_sfcgal(&self) -> Result<SFCGeometry>;
}

/// Convert object to a `CoordSeq` :
/// the (nested) list(s) of coordinates (using tuples of 2 or 3 members)
/// for this geometry.
pub trait ToCoordinates<T> {
    fn to_coordinates(&self) -> Result<CoordSeq<T>>;
}


/// Display SFCGAL version information.
pub fn version() -> String {
    utils::_string(unsafe { sfcgal_version() })
}


#[cfg(test)]
mod tests {
    use super::version;
    #[test]
    fn display_version() {
        assert!(version().contains("1.3."));
    }
}
