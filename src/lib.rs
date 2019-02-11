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
pub use coords::CoordSeq;


/// Convert object to a [`SFCGeometry`] (implemented on [`CoordSeq`] and [geo-types](https://docs.rs/geo-types/) geometries)
///
/// [`SFCGeometry`]: struct.SFCGeometry.html
/// [`CoordSeq`]: enum.CoordSeq.html
pub trait ToSFCGAL {
    fn to_sfcgal(&self) -> Result<SFCGeometry>;
}

/// Convert object to a [`CoordSeq`] holding coordinates and informations about geometry type.
///
/// [`CoordSeq`]: enum.CoordSeq.html
pub trait ToCoordinates {
    // type Item;
    fn to_coordinates<T>(&self) -> Result<CoordSeq<T>> where T: coords::FromSFCGALGeom;
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
