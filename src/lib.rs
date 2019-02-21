//! Rust bindings providing a high-level API to [`SFCGAL`](http://oslandia.github.io/SFCGAL/)
//! library and conversion to / from other geometry crates from Rust ecosystem.
//! Based on the [sfcgal-sys](https://github.com/mthh/sfcgal-rs) crate exposing low-level bindings.
//!
//! Allows notably reading from / writing to WKT as well as interoperability
//! with [geojson]() and [geo-types]() crates.
//! It also offers an API
//! to manipulate SFCGAL geometries from/to coordinates (represented as tuples of 2 or 3 positions).
//!
//! #### Example
//! ```rust
//! # extern crate failure;
//! # fn fun() -> Result<(), failure::Error> {
//! extern crate sfcgal;
//! use sfcgal::{SFCGeometry, CoordSeq, ToGeoJSON, ToSFCGAL, Point2d, Point3d};
//!
//! // Create SFCGAL geometries from coordinates..
//! let coords_linestring = vec![(-0.5, -0.5, 2.5), (0., 0., 4.0)];
//! let coords_polygon = vec![
//!     vec![(-1., -1., 3.0), (1., -1., 3.0), (1., 1., 3.0), (-1., 1., 3.0), (-1., -1., 3.0)], // Exterior ring
//!     vec![(0.1, 0.1, 3.0), (0.1, 0.9, 3.0), (0.9, 0.9, 3.0), (0.9, 0.1, 3.0), (0.1, 0.1, 3.0)], // 1 interior ring
//! ];
//!
//! // .. by using the CoordSeq enum variants to match the wanted SFCGAL geometry type:
//! let line_3d = CoordSeq::Linestring(coords_linestring).to_sfcgal()?;
//! let polygon_3d = CoordSeq::Polygon(coords_polygon).to_sfcgal()?;
//!
//! // ..
//! assert!(line_3d.intersects_3d(&polygon_3d)?);
//! let intersection = line_3d.intersection_3d(&polygon_3d)?;
//!
//! // Get the geojson representation of the geometry with 3-member coordinates :
//! let geom = intersection.to_geojson::<Point3d>()?;
//! // Or the wkt representation with a floating precision of 1 :
//! let wkt = intersection.to_wkt_decim(1)?;
//! # Ok(())
//! # }
//! # fn main() { fun().unwrap(); }
//! ```

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
extern crate geojson;

use sfcgal_sys::sfcgal_version;

mod conversion;
mod errors;
mod utils;
mod geometry;
pub use errors::Result;
pub use geometry::{SFCGeometry, GeomType};
pub use conversion::{CoordSeq, TryInto, AsX3d, X3dString, ToGeoJSON, FromGeoJSON};

/// Type alias for manipulating 2d coordinates, reprensented as (x, y).
pub type Point2d = (f64, f64);

/// Type alias for manipulating 3d coordinates, reprensented as (x, y, z).
pub type Point3d = (f64, f64, f64);

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
    fn to_coordinates<T>(&self) -> Result<CoordSeq<T>> where T: conversion::CoordType + conversion::FromSFCGALGeom;
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
