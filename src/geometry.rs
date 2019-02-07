use sfcgal_sys::{
    initialize,
    sfcgal_io_read_wkt, sfcgal_geometry_as_text, sfcgal_geometry_as_text_decim,
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
    sfcgal_geometry_approximate_medial_axis, sfcgal_geometry_offset_polygon,
};
use std::ffi::{CStr, CString};
use std::ptr::NonNull;
use num_traits::FromPrimitive;
use crate::Result;
use crate::errors::get_last_error;
use crate::utils::{check_predicate, check_computed_value};


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
pub struct SFCGeometry {
    pub(crate) c_geom: NonNull<sfcgal_geometry_t>,
    pub(crate) owned: bool,
}

impl Drop for SFCGeometry {
    fn drop(&mut self) {
        if self.owned == true {
            unsafe { sfcgal_geometry_delete(self.c_geom.as_mut()) }
        }
    }
}

impl Clone for SFCGeometry {
    fn clone(&self) -> SFCGeometry {
        SFCGeometry {
            c_geom: NonNull::new(unsafe { sfcgal_geometry_clone(self.c_geom.as_ref()) }).unwrap(),
            owned: true,
        }
    }
}

impl SFCGeometry {
    pub fn new(wkt: &str) -> Result<SFCGeometry> {
        initialize();
        let c_str = CString::new(wkt)?;
        let obj = unsafe { sfcgal_io_read_wkt(c_str.as_ptr(), wkt.len()) };
        unsafe {
            SFCGeometry::new_from_raw(obj, true)
        }
    }

    pub(crate) unsafe fn new_from_raw(g: *mut sfcgal_geometry_t, owned: bool) -> Result<SFCGeometry> {
        if g.is_null() {
            return Err(
                format_err!(
                 "Reading WKT failed with the following error: {}",
                 get_last_error()
                )
            );
        }
        Ok(
            SFCGeometry {
                c_geom: NonNull::new(g).ok_or(format_err!("Impossible to build the geometry from a Null pointer"))?,
                owned: owned,
            }
        )
    }

    pub fn to_wkt(&self) -> Result<String> {
        let mut ptr: *mut i8 = unsafe { std::mem::uninitialized() };
        let mut length: usize = 0;
        unsafe { sfcgal_geometry_as_text(self.c_geom.as_ref(), &mut ptr, &mut length) };
        let c_str = unsafe { CStr::from_ptr(ptr) };
        Ok(c_str.to_str()?.to_string())
    }

    pub fn to_wkt_decim(&self, nb_decim: i32) -> Result<String> {
        let mut ptr: *mut i8 = unsafe { std::mem::uninitialized() };
        let mut length: usize = 0;
        unsafe { sfcgal_geometry_as_text_decim(self.c_geom.as_ref(), nb_decim, &mut ptr, &mut length) };
        let c_str = unsafe { CStr::from_ptr(ptr) };
        Ok(c_str.to_str()?.to_string())
    }

    pub fn is_empty(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_empty(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_valid(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_valid(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_planar(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_planar(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn is_3d(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_3d(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn _type(&self) -> Result<GeomType> {
        let type_geom = unsafe { sfcgal_geometry_type_id(self.c_geom.as_ptr()) };
        GeomType::from_u32(type_geom)
            .ok_or(format_err!("Unknown geometry type (val={})", type_geom))
   }

    pub fn distance(&self, other: &SFCGeometry) -> Result<f64> {
        let distance = unsafe { sfcgal_geometry_distance(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_computed_value(distance)
    }

    pub fn distance_3d(&self, other: &SFCGeometry) -> Result<f64> {
        let distance = unsafe { sfcgal_geometry_distance_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_computed_value(distance)
    }

    pub fn area(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area(self.c_geom.as_ptr()) };
        check_computed_value(area)
    }

    pub fn area_3d(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area_3d(self.c_geom.as_ptr()) };
        check_computed_value(area)
    }

    pub fn volume(&self) -> Result<f64> {
        let volume = unsafe { sfcgal_geometry_volume(self.c_geom.as_ptr()) };
        check_computed_value(volume)
    }

    pub fn intersects(&self, other: &SFCGeometry) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_intersects(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn intersects_3d(&self, other: &SFCGeometry) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_intersects_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    pub fn intersection(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_intersection(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn intersection_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_intersection_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn difference(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_difference(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn difference_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_difference_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn union(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_union(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn union_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_union_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn straight_skeleton(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_straight_skeleton(self.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn straight_skeleton_distance_in_m(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_straight_skeleton(self.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn approximate_medial_axis(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_approximate_medial_axis(self.c_geom.as_ptr()) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }

    pub fn offset_polygon(&self, radius: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_offset_polygon(self.c_geom.as_ptr(), radius) };
        unsafe {
            SFCGeometry::new_from_raw(result, false)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn geometry_creation_from_wkt() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
    }

    #[test]
    fn geometry_creation_from_wkt_non_valid() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 0.0))");
        assert!(geom.is_ok());
        let geom = geom.unwrap();
        assert!(geom.is_valid().unwrap(), true);
        let geom1 = SFCGeometry::new("POINT(1.0 1.0)").unwrap();
        assert_eq!(geom.intersects(&geom1).unwrap(), true);
    }

    #[test]
    fn geometry_writing_to_wkt() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
        let wkt = geom.unwrap().to_wkt();
        assert!(wkt.is_ok());
        assert_eq!(wkt.unwrap(), String::from("POINT(1/1 1/1)"));
    }

    #[test]
    fn geometry_writing_to_wkt_with_decimals() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
        let wkt = geom.unwrap().to_wkt_decim(1);
        assert!(wkt.is_ok());
        assert_eq!(wkt.unwrap(), String::from("POINT(1.0 1.0)"));
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
