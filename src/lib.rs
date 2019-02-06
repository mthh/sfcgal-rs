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

use std::ffi::CString;
use std::ptr::NonNull;
use sfcgal_sys::{
    initialize,
    sfcgal_io_read_wkt, sfcgal_geometry_as_text,
    sfcgal_geometry_t, sfcgal_geometry_clone, sfcgal_geometry_delete, sfcgal_geometry_type_id,
    sfcgal_geometry_is_empty, sfcgal_geometry_is_3d, sfcgal_geometry_is_planar,
    sfcgal_geometry_is_valid, sfcgal_geometry_is_valid_detail,
    sfcgal_geometry_area, sfcgal_geometry_area_3d, sfcgal_geometry_volume,
    sfcgal_geometry_distance, sfcgal_geometry_distance_3d,
    sfcgal_geometry_intersects, sfcgal_geometry_intersects_3d,
    sfcgal_geometry_intersection, sfcgal_geometry_intersection_3d,
    sfcgal_geometry_difference, sfcgal_geometry_difference_3d,
    sfcgal_geometry_union, sfcgal_geometry_union_3d,
    sfcgal_geometry_straight_skeleton, sfcgal_geometry_straight_skeleton_distance_in_m,
    sfcgal_prepared_geometry_delete,
    };
use num_traits::FromPrimitive;

mod geo;
mod errors;
mod utils;
pub use errors::Result;

#[repr(C)]
#[derive(PartialEq, Eq, Debug, Primitive)]
pub enum GeomType {
  Point = 1,
  Linestring = 2,
  Polygon = 3,
  Multipoint = 4,
  Multilinestring = 5,
  Multipolygon = 6,
  Geometrycollection = 7,
  Polyhedralsurface = 15,
  Triangulatedsurface = 16,
  Triangle = 100,
  Solid = 101,
  Multisolid = 102
}

#[repr(C)]
pub struct SFCGeometry(NonNull<sfcgal_geometry_t>);

impl Drop for SFCGeometry {
    fn drop(&mut self) {
        unsafe { sfcgal_geometry_delete(self.0.as_mut()) }
    }
}

impl Clone for SFCGeometry {
    fn clone(&self) -> SFCGeometry {
        SFCGeometry(NonNull::new(unsafe { sfcgal_geometry_clone(self.0.as_ref()) }).unwrap())
    }
}

/// Convert object to a SFCGAL geometry.
pub trait ToSfcgal {
    fn to_sfcgal(&self) -> Result<SFCGeometry>;
}

impl SFCGeometry {
    pub fn new(wkt: &str) -> Result<SFCGeometry> {
        initialize();
        let c_str = CString::new(wkt)?;
        let obj = unsafe { sfcgal_io_read_wkt(c_str.as_ptr(), wkt.len()) };
        unsafe {
            SFCGeometry::new_from_raw(obj)
        }
    }

    unsafe fn new_from_raw(g: *mut sfcgal_geometry_t) -> Result<SFCGeometry> {
        if g.is_null() {
            return Err(format_err!(
                 "Reading WKT failed with the following error: {}",
                 errors::get_last_error()
            ));
        }
        NonNull::new(g)
            .ok_or(format_err!("Impossible to build the geometry from a nullptr"))
            .map(SFCGeometry)
    }

    pub fn to_wkt(&self) -> Result<String> {
        // let mut c_wkt = std::ptr::null_mut();
        // let mut length: *mut usize = std::ptr::null_mut::<usize>();
        // unsafe { sfcgal_geometry_as_text(self.0.as_ptr(), c_wkt as *mut *mut i8, length) };
        // let wkt = utils::_string(c_wkt);
        // Ok(wkt)
        let c_string = CString::new("").expect("CString::new failed");
        let length: *mut usize = std::ptr::null_mut::<usize>();
        let raw = c_string.into_raw();
        unsafe {
            sfcgal_geometry_as_text(self.0.as_ptr(), raw as *mut *mut i8, length)
        };
        let res = unsafe { CString::from_raw(raw).into_string()
            .expect("Error while creating the WKT representation of the given geometry") };
        Ok(res)
    }

    pub fn is_empty(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_empty(self.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_valid(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_valid(self.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_planar(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_planar(self.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_3d(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_3d(self.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn _type(&self) -> Result<GeomType> {
        let type_geom = unsafe { sfcgal_geometry_type_id(self.0.as_ptr()) };
        GeomType::from_u32(type_geom)
            .ok_or(format_err!("Unknown geometry type (val={})", type_geom))
   }

    pub fn distance(&self, other: &SFCGeometry) -> Result<f64> {
        let distance = unsafe { sfcgal_geometry_distance(self.0.as_ptr(), other.0.as_ptr()) };
        check_computed_value(distance)
    }

    pub fn distance_3d(&self, other: &SFCGeometry) -> Result<f64> {
        let distance = unsafe { sfcgal_geometry_distance_3d(self.0.as_ptr(), other.0.as_ptr()) };
        check_computed_value(distance)
    }

    pub fn area(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area(self.0.as_ptr()) };
        check_computed_value(area)
    }

    pub fn area_3d(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area_3d(self.0.as_ptr()) };
        check_computed_value(area)
    }

    pub fn volume(&self) -> Result<f64> {
        let volume = unsafe { sfcgal_geometry_area_3d(self.0.as_ptr()) };
        check_computed_value(volume)
    }

    pub fn intersects(&self, other: &SFCGeometry) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_intersects(self.0.as_ptr(), other.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn intersects_3d(&self, other: &SFCGeometry) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_intersects_3d(self.0.as_ptr(), other.0.as_ptr()) };
        check_predicate(rv)
    }

    pub fn intersection(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_intersection(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn intersection_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_intersection_3d(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn difference(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_difference(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn difference_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_difference_3d(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn union(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_union(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn union_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_union_3d(self.0.as_ptr(), other.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }

    pub fn straight_skeleton(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_straight_skeleton(self.0.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result)
        }
    }
}

fn check_predicate(val: i32) -> Result<bool> {
    match val {
        1 => Ok(true),
        0 => Ok(false),
        _ => Err(format_err!("SFCGAL error: {}", errors::get_last_error()))
    }
}

fn check_computed_value(val: f64) -> Result<f64> {
    match val {
        -1.0 => Err(format_err!("SFCGAL error: {}", errors::get_last_error())),
        _ => Ok(val)
    }
}

#[repr(C)]
pub struct SFCPreparedGeometry(NonNull<sfcgal_geometry_t>);

impl Drop for SFCPreparedGeometry {
    fn drop(&mut self) {
        unsafe { sfcgal_prepared_geometry_delete(self.0.as_mut()) }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn geometry_creation() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
    }

    #[test]
    fn geometry_creation_failed_with_error() {
        let geom = SFCGeometry::new("POINT(1, 1)");
        assert!(geom.is_err());
        assert_eq!(
            geom.err().unwrap().to_string(),
            "Reading WKT failed with the following error: WKT parse error, Coordinate dimension < 2 (, 1))",
        )
    }

    #[test]
    fn geometry_distance_to_other() {
        let pt1 = SFCGeometry::new("POINT(1.0 1.0)").unwrap();
        let pt2 = SFCGeometry::new("POINT(10.0 1.0)").unwrap();
        let distance = pt1.distance(&pt2).unwrap();
        assert_eq!(distance, 9.0);
    }

    #[test]
    fn geometry_distance_3d_to_other() {
        let pt1 = SFCGeometry::new("POINT(1.0 1.0 2.0)").unwrap();
        let pt2 = SFCGeometry::new("POINT(10.0 1.0 2.0)").unwrap();
        let distance = pt1.distance_3d(&pt2).unwrap();
        assert_eq!(distance, 9.0);
    }

    #[test]
    fn geometry_area() {
        let polygon = SFCGeometry::new("POLYGON((1 1, 3 1, 4 4, 1 3, 1 1))").unwrap();
        assert_eq!(polygon.area().unwrap(), 6.0);
    }

    #[test]
    fn geometry_volume() {
        let surface = SFCGeometry::new("POLYHEDRALSURFACE Z (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)),((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),((1 1 1, 1 0 1, 0 0 1, 0 1 1, 1 1 1)),((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))").unwrap();
        assert_eq!(surface.volume().is_err(), true);
    }

    #[test]
    fn geometry_predicates() {
        let pt = SFCGeometry::new("POINT(1.0 1.0)").unwrap();
        assert_eq!(pt.is_valid().unwrap(), true);
        assert_eq!(pt.is_3d().unwrap(), false);
        assert_eq!(pt.is_empty().unwrap(), false);
        assert_eq!(
            pt.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let linestring_3d = SFCGeometry::new("LINESTRING(10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();
        assert_eq!(linestring_3d.is_valid().unwrap(), true);
        assert_eq!(linestring_3d.is_3d().unwrap(), true);
        assert_eq!(linestring_3d.is_empty().unwrap(), false);
        assert_eq!(
            linestring_3d.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let empty_geom = SFCGeometry::new("LINESTRING EMPTY").unwrap();
        assert_eq!(empty_geom.is_valid().unwrap(), true);
        assert_eq!(empty_geom.is_3d().unwrap(), false);
        assert_eq!(empty_geom.is_empty().unwrap(), true);
        assert_eq!(
            linestring_3d.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let polyg = SFCGeometry::new("POLYGON((1 1, 3 1, 4 4, 1 3, 1 1))").unwrap();
        assert_eq!(polyg.is_valid().unwrap(), true);
        assert_eq!(polyg.is_3d().unwrap(), false);
        assert_eq!(polyg.is_empty().unwrap(), false);
        assert_eq!(polyg.is_planar().unwrap(), true);

        assert_eq!(pt.intersects(&polyg).unwrap(), true);
        assert_eq!(pt.intersects_3d(&linestring_3d).unwrap(), false);
    }
}
