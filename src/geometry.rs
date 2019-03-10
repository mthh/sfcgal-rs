use crate::conversion::{CoordSeq, CoordType, ToSFCGALGeom};
use crate::errors::get_last_error;
use crate::utils::{_c_string_with_size, _string, check_computed_value, check_predicate};
use crate::{Result, ToSFCGAL};
use num_traits::FromPrimitive;
use sfcgal_sys::{
    initialize, sfcgal_geometry_approximate_medial_axis, sfcgal_geometry_area,
    sfcgal_geometry_area_3d, sfcgal_geometry_as_text, sfcgal_geometry_as_text_decim,
    sfcgal_geometry_clone, sfcgal_geometry_convexhull, sfcgal_geometry_convexhull_3d,
    sfcgal_geometry_delete, sfcgal_geometry_difference, sfcgal_geometry_difference_3d,
    sfcgal_geometry_distance, sfcgal_geometry_distance_3d, sfcgal_geometry_extrude,
    sfcgal_geometry_intersection, sfcgal_geometry_intersection_3d, sfcgal_geometry_intersects,
    sfcgal_geometry_intersects_3d, sfcgal_geometry_is_3d, sfcgal_geometry_is_empty,
    sfcgal_geometry_is_planar, sfcgal_geometry_is_valid, sfcgal_geometry_is_valid_detail,
    sfcgal_geometry_offset_polygon, sfcgal_geometry_orientation, sfcgal_geometry_straight_skeleton,
    sfcgal_geometry_straight_skeleton_distance_in_m, sfcgal_geometry_t, sfcgal_geometry_tesselate,
    sfcgal_geometry_triangulate_2dz, sfcgal_geometry_type_id, sfcgal_geometry_union,
    sfcgal_geometry_union_3d, sfcgal_geometry_volume, sfcgal_io_read_wkt,
};
use std::ffi::CString;
use std::ptr::NonNull;

/// SFCGAL Geometry types.
///
/// Indicates the type of shape represented by a `SFCGeometry`.
/// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga1afcf1fad6c2daeca001481b125b84c6))
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
    Multisolid = 102,
}

/// Represents the orientation of a `SFCGeometry`.
#[derive(PartialEq, Eq, Debug, Primitive)]
pub enum Orientation {
    CounterClockWise = -1isize,
    ClockWise = 1isize,
    Undetermined = 0isize,
}

/// Object representing a SFCGAL Geometry.
///
/// Most of the operations allowed by SFCGAL C API are wrapped,
/// except those modifying the geometry in-place (such as adding a new
/// point to a linestring for example) and those retrieving a specific part
/// of a geometry (such as getting the 2nd interior ring of some polygon as a linestring).
/// However, this can easily be done by yourself by converting them from/to coordinates
/// with the `new_from_coordinates` and `to_coordinates` methods.
///
/// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#gadd6d3ea5a71a957581248791624fad58))
#[repr(C)]
pub struct SFCGeometry {
    pub(crate) c_geom: NonNull<sfcgal_geometry_t>,
    pub(crate) owned: bool,
}

impl Drop for SFCGeometry {
    fn drop(&mut self) {
        if self.owned {
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

impl std::fmt::Debug for SFCGeometry {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.to_wkt_decim(8).unwrap())
    }
}

impl SFCGeometry {
    /// Create a geometry by parsing a [WKT](https://en.wikipedia.org/wiki/Well-known_text) string.
    pub fn new(wkt: &str) -> Result<SFCGeometry> {
        initialize();
        let c_str = CString::new(wkt)?;
        let obj = unsafe { sfcgal_io_read_wkt(c_str.as_ptr(), wkt.len()) };
        unsafe { SFCGeometry::new_from_raw(obj, true) }
    }

    pub(crate) unsafe fn new_from_raw(
        g: *mut sfcgal_geometry_t,
        owned: bool,
    ) -> Result<SFCGeometry> {
        Ok(SFCGeometry {
            owned,
            c_geom: NonNull::new(g).ok_or_else(|| {
                format_err!(
                    "Obtained null pointer when creating geometry: {}",
                    get_last_error()
                )
            })?,
        })
    }

    pub fn new_from_coordinates<T>(coords: &CoordSeq<T>) -> Result<SFCGeometry>
    where
        T: ToSFCGALGeom + CoordType,
    {
        coords.to_sfcgal()
    }

    /// Returns a WKT representation of the given `SFCGeometry` using CGAL exact integer fractions as coordinate values.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga3bc1954e3c034b60f0faff5e8227c398))
    pub fn to_wkt(&self) -> Result<String> {
        let mut ptr: *mut i8 = unsafe { std::mem::uninitialized() };
        let mut length: usize = 0;
        unsafe { sfcgal_geometry_as_text(self.c_geom.as_ref(), &mut ptr, &mut length) };
        Ok(_c_string_with_size(ptr, length))
    }

    /// Returns a WKT representation of the given `SFCGeometry` using floating point coordinate values with
    /// the desired number of decimals.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#gaaf23f2c95fd48810beb37d07a9652253))
    pub fn to_wkt_decim(&self, nb_decim: i32) -> Result<String> {
        let mut ptr: *mut i8 = unsafe { std::mem::uninitialized() };
        let mut length: usize = 0;
        unsafe {
            sfcgal_geometry_as_text_decim(self.c_geom.as_ref(), nb_decim, &mut ptr, &mut length)
        };
        Ok(_c_string_with_size(ptr, length))
    }

    /// Test if the given `SFCGeometry` is empty or not.
    pub fn is_empty(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_empty(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Test if the given `SFCGeometry` is valid or not.
    pub fn is_valid(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_valid(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Returns reason for the invalidity or None in case of validity.
    pub fn validity_detail(&self) -> Result<Option<String>> {
        let mut ptr: *mut i8 = unsafe { std::mem::uninitialized() };
        let rv = unsafe {
            sfcgal_geometry_is_valid_detail(
                self.c_geom.as_ptr(),
                &mut ptr,
                std::ptr::null::<sfcgal_geometry_t>() as *mut *mut sfcgal_geometry_t,
            )
        };
        match rv {
            1 => Ok(None),
            0 => Ok(Some(_string(ptr))),
            _ => Err(format_err!("SFCGAL error: {}", get_last_error())),
        }
    }

    /// Test if the given `SFCGeometry` is planar or not.
    pub fn is_planar(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_planar(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Test if the given `SFCGeometry` is a 3d geometry or not.
    pub fn is_3d(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_3d(self.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Returns the SFCGAL type of the given `SFCGeometry`.
    pub fn _type(&self) -> Result<GeomType> {
        let type_geom = unsafe { sfcgal_geometry_type_id(self.c_geom.as_ptr()) };
        GeomType::from_u32(type_geom)
            .ok_or_else(|| format_err!("Unknown geometry type (val={})", type_geom))
    }

    /// Computes the distance to an other `SFCGeometry`.
    pub fn distance(&self, other: &SFCGeometry) -> Result<f64> {
        let distance =
            unsafe { sfcgal_geometry_distance(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_computed_value(distance)
    }

    /// Computes the 3d distance to an other `SFCGeometry`.
    pub fn distance_3d(&self, other: &SFCGeometry) -> Result<f64> {
        let distance =
            unsafe { sfcgal_geometry_distance_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_computed_value(distance)
    }

    /// Computes the area of the given `SFCGeometry`.
    pub fn area(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area(self.c_geom.as_ptr()) };
        check_computed_value(area)
    }

    /// Computes the 3d area of the given `SFCGeometry`.
    pub fn area_3d(&self) -> Result<f64> {
        let area = unsafe { sfcgal_geometry_area_3d(self.c_geom.as_ptr()) };
        check_computed_value(area)
    }

    /// Computes the volume of the given `SFCGeometry` (must be a volume).
    pub fn volume(&self) -> Result<f64> {
        let volume = unsafe { sfcgal_geometry_volume(self.c_geom.as_ptr()) };
        check_computed_value(volume)
    }

    /// Computes the orientation of the given `SFCGeometry` (must be a Polygon)
    pub fn orientation(&self) -> Result<Orientation> {
        let orientation = unsafe { sfcgal_geometry_orientation(self.c_geom.as_ptr()) };
        Orientation::from_i32(orientation)
            .ok_or_else(|| format_err!("Error while retrieving orientation (val={})", orientation))
    }

    /// Test the intersection with an other `SFCGeometry`.
    pub fn intersects(&self, other: &SFCGeometry) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_intersects(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Test the 3d intersection with an other `SFCGeometry`.
    pub fn intersects_3d(&self, other: &SFCGeometry) -> Result<bool> {
        let rv =
            unsafe { sfcgal_geometry_intersects_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        check_predicate(rv)
    }

    /// Returns the intersection of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn intersection(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_intersection(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d intersection of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn intersection_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_intersection_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the difference of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn difference(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_difference(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d difference of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn difference_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_difference_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the union of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn union(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_union(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d union of the given `SFCGeometry` to an other `SFCGeometry`.
    pub fn union_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_union_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the straight skeleton of the given `SFCGeometry`.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#gaefaa76b61d66e2ad11d902e6b5a13635))
    pub fn straight_skeleton(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_straight_skeleton(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }
    /// Returns the straight skeleton of the given `SFCGeometry` with the distance to the border as M coordinate.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga972ea9e378eb2dc99c00b6ad57d05e88))
    pub fn straight_skeleton_distance_in_m(&self) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_straight_skeleton_distance_in_m(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the approximate medial axis for the given `SFCGeometry` Polygon.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga16a9b4b1211843f8444284b1fefebc46))
    pub fn approximate_medial_axis(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_approximate_medial_axis(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the offset polygon of the given `SFCGeometry`.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga9766f54ebede43a9b71fccf1524a1054))
    pub fn offset_polygon(&self, radius: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_offset_polygon(self.c_geom.as_ptr(), radius) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the extrusion of the given `SFCGeometry` (not supported on Solid and Multisolid).
    /// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga277d01bd9978e13644baa1755f1cd3e0)
    pub fn extrude(&self, ex: f64, ey: f64, ez: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_extrude(self.c_geom.as_ptr(), ex, ey, ez) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns a tesselation of the given `SFCGeometry`.
    /// ([C API reference](http://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga570ce6214f305ed35ebbec62d366b588))
    pub fn tesselate(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_tesselate(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns a triangulation of the given `SFCGeometry`.
    /// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#gae382792f387654a9adb2e2c38735e08d))
    pub fn triangulate_2dz(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_triangulate_2dz(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the convex hull of the given `SFCGeometry`.
    /// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#ga9027b5654cbacf6c2106d70b129d3a23))
    pub fn convexhull(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_convexhull(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d convex hull of the given `SFCGeometry`.
    /// ([C API reference](https://oslandia.github.io/SFCGAL/doxygen/group__capi.html#gacf01a9097f2059afaad871658b4b5a6f))
    pub fn convexhull_3d(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_convexhull_3d(self.c_geom.as_ptr()) };
        unsafe { SFCGeometry::new_from_raw(result, true) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creation_point_from_wkt() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
    }

    #[test]
    fn creation_polygon_from_wkt() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 0.0))");
        assert!(geom.is_ok());
        let geom = geom.unwrap();
        assert!(geom.is_valid().unwrap(), true);
        let geom1 = SFCGeometry::new("POINT(1.0 1.0)").unwrap();
        assert_eq!(geom.intersects(&geom1).unwrap(), true);
    }

    #[test]
    fn writing_to_wkt() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
        let wkt = geom.unwrap().to_wkt();
        assert!(wkt.is_ok());
        assert_eq!(wkt.unwrap(), String::from("POINT(1/1 1/1)"));
    }

    #[test]
    fn writing_to_wkt_with_decimals() {
        let geom = SFCGeometry::new("POINT(1.0 1.0)");
        assert!(geom.is_ok());
        let wkt = geom.unwrap().to_wkt_decim(1);
        assert!(wkt.is_ok());
        assert_eq!(wkt.unwrap(), String::from("POINT(1.0 1.0)"));
    }

    #[test]
    fn creation_failed_with_error_message() {
        let geom = SFCGeometry::new("POINT(1, 1)");
        assert!(geom.is_err());
        assert_eq!(
            geom.err().unwrap().to_string(),
            "Obtained null pointer when creating geometry: WKT parse error, Coordinate dimension < 2 (, 1))",
        )
    }

    #[test]
    fn distance_to_other() {
        let pt1 = SFCGeometry::new("POINT(1.0 1.0)").unwrap();
        let pt2 = SFCGeometry::new("POINT(10.0 1.0)").unwrap();
        let distance = pt1.distance(&pt2).unwrap();
        assert_eq!(distance, 9.0);
    }

    #[test]
    fn distance_3d_to_other() {
        let pt1 = SFCGeometry::new("POINT(1.0 1.0 2.0)").unwrap();
        let pt2 = SFCGeometry::new("POINT(10.0 1.0 2.0)").unwrap();
        let distance = pt1.distance_3d(&pt2).unwrap();
        assert_eq!(distance, 9.0);
    }

    #[test]
    fn area() {
        let polygon = SFCGeometry::new("POLYGON((1 1, 3 1, 4 4, 1 3, 1 1))").unwrap();
        assert_eq!(polygon.area().unwrap(), 6.0);
    }

    #[test]
    fn area_3d() {
        let polygon = SFCGeometry::new("POLYGON((1 1 1, 3 1 1, 4 4 1, 1 3 1, 1 1 1))").unwrap();
        assert_ulps_eq!(polygon.area_3d().unwrap(), 6.0);
    }

    #[test]
    fn volume() {
        let cube = SFCGeometry::new(
            "SOLID((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
             ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
             ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
             ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
             ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
             ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))",
        )
        .unwrap();
        assert_eq!(cube.volume().unwrap(), 1.);
    }

    #[test]
    fn volume_on_not_volume_geometry() {
        let surface = SFCGeometry::new(
            "POLYHEDRALSURFACE Z \
             (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
             ((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)),\
             ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
             ((1 1 1, 1 0 1, 0 0 1, 0 1 1, 1 1 1)),\
             ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
             ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))",
        )
        .unwrap();
        assert_eq!(surface.volume().is_err(), true);
    }

    #[test]
    fn predicates() {
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

    #[test]
    fn validity_detail_on_valid_geom() {
        let line = SFCGeometry::new("LINESTRING(10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();
        assert_eq!(line.is_valid().unwrap(), true);
        assert_eq!(line.validity_detail().unwrap(), None);
    }

    #[test]
    fn validity_detail_on_invalid_geom() {
        let surface = SFCGeometry::new(
            "POLYHEDRALSURFACE Z \
             (((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
             ((0 0 0, 0 1 0, 0 1 1, 0 0 1, 0 0 0)),\
             ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
             ((1 1 1, 1 0 1, 0 0 1, 0 1 1, 1 1 1)),\
             ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
             ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1)))",
        )
        .unwrap();
        assert_eq!(surface.is_valid().unwrap(), false);
        assert_eq!(
            surface.validity_detail().unwrap(),
            Some(String::from("inconsistent orientation of PolyhedralSurface detected at edge 3 (4-7) of polygon 5")),
        );
    }

    #[test]
    fn straight_skeleton() {
        let geom = SFCGeometry::new("POLYGON((0 0,1 0,1 1,0 1,0 0))").unwrap();
        let result = geom.straight_skeleton().unwrap();
        let wkt = result.to_wkt_decim(1).unwrap();
        assert_eq!(
            wkt,
            "MULTILINESTRING((-0.0 -0.0,0.5 0.5),(1.0 -0.0,0.5 0.5),(1.0 1.0,0.5 0.5),(-0.0 1.0,0.5 0.5))",
        );
    }

    #[test]
    fn straight_skeleton_distance_in_m() {
        let geom = SFCGeometry::new("POLYGON((0 0,1 0,1 1,0 1,0 0))").unwrap();
        let result = geom.straight_skeleton_distance_in_m().unwrap();
        let wkt = result.to_wkt_decim(1).unwrap();
        assert_eq!(
            wkt,
            "MULTILINESTRING M((-0.0 -0.0 0.0,0.5 0.5 0.5),(1.0 -0.0 0.0,0.5 0.5 0.5),(1.0 1.0 0.0,0.5 0.5 0.5),(-0.0 1.0 0.0,0.5 0.5 0.5))",
        );
    }

    #[test]
    fn tesselate() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();
        let result = geom.tesselate().unwrap();
        let output_wkt = result.to_wkt_decim(1).unwrap();
        assert_eq!(
            output_wkt,
            "TIN(((0.0 1.0,1.0 0.0,1.0 1.0,0.0 1.0)),((0.0 1.0,0.0 0.0,1.0 0.0,0.0 1.0)))",
        );
    }

    #[test]
    fn offset_polygon() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();
        let buff = geom.offset_polygon(1.).unwrap();
        assert!(buff.is_valid().unwrap());
        assert!(!buff.is_empty().unwrap());
    }

    #[test]
    fn extrude_polygon() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();
        let extr = geom.extrude(0., 0., 1.).unwrap();
        assert!(extr.is_valid().unwrap());
        assert!(!extr.is_empty().unwrap());
        assert_eq!(extr._type().unwrap(), GeomType::Solid);
    }

    #[test]
    fn tesselate_invariant_geom() {
        let input_wkt = String::from("POINT(1.0 1.0)");
        let pt = SFCGeometry::new(&input_wkt).unwrap();
        let result = pt.tesselate().unwrap();
        let output_wkt = result.to_wkt_decim(1).unwrap();
        assert_eq!(input_wkt, output_wkt);
    }

    #[test]
    fn difference_3d() {
        let cube1 = SFCGeometry::new(
            "
            SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
            ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
            ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
            ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
            ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
            ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))))",
        )
        .unwrap();
        let cube2 = SFCGeometry::new(
            "
            SOLID((((0 0 0.5, 0 1 0.5, 1 1 0.5, 1 0 0.5, 0 0 0.5)),\
            ((0 0 0.5, 0 0 1, 0 1 1, 0 1 0.5, 0 0 0.5)),\
            ((0 0 0.5, 1 0 0.5, 1 0 1, 0 0 1, 0 0 0.5)),\
            ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
            ((1 1 1, 1 0 1, 1 0 0.5, 1 1 0.5, 1 1 1)),\
            ((1 1 1, 1 1 0.5, 0 1 0.5, 0 1 1, 1 1 1))))",
        )
        .unwrap();
        let diff = cube1.difference_3d(&cube2).unwrap();
        assert_eq!(diff.is_valid().unwrap(), true);
        assert_ulps_eq!(diff.volume().unwrap(), 0.5);
    }

    #[test]
    fn intersection_3d() {
        let cube1 = SFCGeometry::new(
            "
            SOLID((((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
            ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
            ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
            ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
            ((1 1 1, 1 0 1, 1 0 0, 1 1 0, 1 1 1)),\
            ((1 1 1, 1 1 0, 0 1 0, 0 1 1, 1 1 1))))",
        )
        .unwrap();
        let cube2 = SFCGeometry::new(
            "
            SOLID((((0 0 0.5, 0 1 0.5, 1 1 0.5, 1 0 0.5, 0 0 0.5)),\
            ((0 0 0.5, 0 0 1, 0 1 1, 0 1 0.5, 0 0 0.5)),\
            ((0 0 0.5, 1 0 0.5, 1 0 1, 0 0 1, 0 0 0.5)),\
            ((1 1 1, 0 1 1, 0 0 1, 1 0 1, 1 1 1)),\
            ((1 1 1, 1 0 1, 1 0 0.5, 1 1 0.5, 1 1 1)),\
            ((1 1 1, 1 1 0.5, 0 1 0.5, 0 1 1, 1 1 1))))",
        )
        .unwrap();
        let diff = cube1.intersection_3d(&cube2).unwrap();
        assert_eq!(diff.is_valid().unwrap(), true);
        assert_ulps_eq!(diff.volume().unwrap(), 0.5);
    }
}
