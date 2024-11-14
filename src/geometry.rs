use std::{
    ffi::{c_void, CString},
    mem::MaybeUninit,
    os::raw::c_char,
    ptr::NonNull,
};

use num_traits::FromPrimitive;
use sfcgal_sys::{
    initialize, sfcgal_alloc_handler_t, sfcgal_error_handler_t, sfcgal_free_handler_t,
    sfcgal_geometry_alpha_shapes, sfcgal_geometry_approximate_medial_axis, sfcgal_geometry_area,
    sfcgal_geometry_area_3d, sfcgal_geometry_as_text, sfcgal_geometry_as_text_decim,
    sfcgal_geometry_clone, sfcgal_geometry_collection_add_geometry,
    sfcgal_geometry_collection_create, sfcgal_geometry_collection_geometry_n,
    sfcgal_geometry_collection_num_geometries, sfcgal_geometry_convexhull,
    sfcgal_geometry_convexhull_3d, sfcgal_geometry_covers, sfcgal_geometry_covers_3d,
    sfcgal_geometry_delete, sfcgal_geometry_difference, sfcgal_geometry_difference_3d,
    sfcgal_geometry_distance, sfcgal_geometry_distance_3d, sfcgal_geometry_extrude,
    sfcgal_geometry_extrude_polygon_straight_skeleton, sfcgal_geometry_extrude_straight_skeleton,
    sfcgal_geometry_intersection, sfcgal_geometry_intersection_3d, sfcgal_geometry_intersects,
    sfcgal_geometry_intersects_3d, sfcgal_geometry_is_3d, sfcgal_geometry_is_empty,
    sfcgal_geometry_is_measured, sfcgal_geometry_is_planar, sfcgal_geometry_is_valid,
    sfcgal_geometry_is_valid_detail, sfcgal_geometry_line_sub_string,
    sfcgal_geometry_minkowski_sum, sfcgal_geometry_offset_polygon,
    sfcgal_geometry_optimal_alpha_shapes, sfcgal_geometry_orientation,
    sfcgal_geometry_straight_skeleton, sfcgal_geometry_straight_skeleton_distance_in_m,
    sfcgal_geometry_t, sfcgal_geometry_tesselate, sfcgal_geometry_triangulate_2dz,
    sfcgal_geometry_type_id, sfcgal_geometry_union, sfcgal_geometry_union_3d,
    sfcgal_geometry_volume, sfcgal_io_read_wkt, sfcgal_multi_linestring_create,
    sfcgal_multi_point_create, sfcgal_multi_polygon_create, sfcgal_prepared_geometry_t, srid_t,
};
use sfcgal_sys::{
    /* sfcgal_solid_set_exterior_shell, */
    sfcgal_approx_convex_partition_2, sfcgal_full_version, sfcgal_geometry_as_hexwkb,
    sfcgal_geometry_as_obj, sfcgal_geometry_as_obj_file, sfcgal_geometry_as_vtk,
    sfcgal_geometry_as_vtk_file, sfcgal_geometry_as_wkb, sfcgal_geometry_buffer3d,
    sfcgal_geometry_force_lhr, sfcgal_geometry_force_rhr, sfcgal_geometry_force_valid,
    sfcgal_geometry_has_validity_flag, sfcgal_geometry_make_solid, sfcgal_geometry_rotate,
    sfcgal_geometry_rotate_2d, sfcgal_geometry_rotate_3d, sfcgal_geometry_rotate_3d_around_center,
    sfcgal_geometry_rotate_x, sfcgal_geometry_rotate_y, sfcgal_geometry_rotate_z,
    sfcgal_geometry_round, sfcgal_geometry_scale, sfcgal_geometry_scale_3d,
    sfcgal_geometry_scale_3d_around_center, sfcgal_geometry_straight_skeleton_partition,
    sfcgal_geometry_translate_2d, sfcgal_geometry_translate_3d, sfcgal_geometry_visibility_point,
    sfcgal_geometry_visibility_segment, sfcgal_greene_approx_convex_partition_2, sfcgal_init,
    sfcgal_io_read_binary_prepared, sfcgal_io_read_ewkt, sfcgal_io_read_wkb,
    sfcgal_io_write_binary_prepared, sfcgal_linestring_add_point, sfcgal_linestring_create,
    sfcgal_linestring_num_points, sfcgal_linestring_point_n, sfcgal_multi_solid_create,
    sfcgal_optimal_convex_partition_2, sfcgal_point_create, sfcgal_point_create_from_xy,
    sfcgal_point_create_from_xym, sfcgal_point_create_from_xyz, sfcgal_point_create_from_xyzm,
    sfcgal_point_m, sfcgal_point_x, sfcgal_point_y, sfcgal_point_z,
    sfcgal_polygon_add_interior_ring, sfcgal_polygon_create,
    sfcgal_polygon_create_from_exterior_ring, sfcgal_polygon_exterior_ring,
    sfcgal_polygon_interior_ring_n, sfcgal_polygon_num_interior_rings,
    sfcgal_polyhedral_surface_add_polygon, sfcgal_polyhedral_surface_create,
    sfcgal_polyhedral_surface_num_polygons, sfcgal_polyhedral_surface_polygon_n,
    sfcgal_prepared_geometry_as_ewkt, sfcgal_prepared_geometry_create,
    sfcgal_prepared_geometry_create_from_geometry, sfcgal_prepared_geometry_delete,
    sfcgal_prepared_geometry_geometry, sfcgal_prepared_geometry_set_geometry,
    sfcgal_prepared_geometry_set_srid, sfcgal_prepared_geometry_srid, sfcgal_set_alloc_handlers,
    sfcgal_set_error_handlers, sfcgal_set_geometry_validation, sfcgal_solid_add_interior_shell,
    sfcgal_solid_create, sfcgal_solid_create_from_exterior_shell, sfcgal_solid_num_shells,
    sfcgal_solid_shell_n, sfcgal_triangle_create, sfcgal_triangle_create_from_points,
    sfcgal_triangle_set_vertex, sfcgal_triangle_set_vertex_from_xy,
    sfcgal_triangle_set_vertex_from_xyz, sfcgal_triangle_vertex,
    sfcgal_triangulated_surface_add_triangle, sfcgal_triangulated_surface_create,
    sfcgal_triangulated_surface_num_triangles, sfcgal_triangulated_surface_triangle_n,
    sfcgal_version, sfcgal_y_monotone_partition_2,
};

use crate::{
    conversion::{CoordSeq, CoordType, ToSFCGALGeom},
    errors::get_last_error,
    utils::{
        _c_string_with_size, _string, check_computed_value, check_nan_value,
        check_null_prepared_geom, check_predicate, get_raw_bytes,
    },
    Result, ToSFCGAL,
};

#[repr(C)]
pub enum BufferType {
    Round,
    CylSphere,
    Flat,
}

/// SFCGAL Geometry types.
///
/// Indicates the type of shape represented by a `SFCGeometry`.
#[repr(C)]
#[derive(PartialEq, Eq, PartialOrd, Ord, Debug, Primitive)]

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
    Triangle = 17,
    Solid = 101,
    Multisolid = 102,
}

impl GeomType {
    fn is_collection_type(&self) -> bool {
        matches!(
            &self,
            GeomType::Multipoint
                | GeomType::Multilinestring
                | GeomType::Multipolygon
                | GeomType::Multisolid
                | GeomType::Geometrycollection
        )
    }
}

/// Represents the orientation of a `SFCGeometry`.
#[derive(PartialEq, Eq, Debug, Primitive)]

pub enum Orientation {
    CounterClockWise = -1isize,
    ClockWise = 1isize,
    Undetermined = 0isize,
}

macro_rules! precondition_match_type {
    ($value : expr, $input_type: expr) => {{
        if $value._type()? != $input_type {
            bail!("Wrong input Geometry type")
        }
    }};
}

macro_rules! precondition_match_type_other {
    ($value : expr, $input_type: expr) => {{
        if $value._type()? != $input_type {
            bail!("Wrong input Geometry type for other")
        }
    }};
}

macro_rules! precondition_match_validity {
    ($value : expr) => {{
        if !$value.is_valid()? {
            bail!("Input geometry is not valid")
        }
    }};
}

macro_rules! precondition_match_validity_other {
    ($value : expr) => {{
        if !$value.is_valid()? {
            bail!("Input geometry is not valid for other")
        }
    }};
}

macro_rules! precondition_match_not_empty {
    ($value : expr) => {{
        if $value.is_empty()? {
            bail!("Input geometry is empty")
        }
    }};
}

macro_rules! precondition_index_in_range {
    ($index: expr, $range: expr) => {{
        if !$range.contains(&$index) {
            bail!("Index not in the expected range")
        }
    }};
}

macro_rules! precondition_index_in_result_value {
    ($index: expr, $range_result: expr) => {{
        match $range_result {
            Ok(max_value) => {
                if !(0..max_value).contains(&$index) {
                    bail!("Index not in the expected range")
                }
            }
            Err(err) => {
                bail!(err)
            }
        }
    }};
}

/// Object representing a SFCGAL Geometry.
///
/// Most of the operations allowed by SFCGAL C API are wrapped,
/// except those modifying the geometry in-place (such as adding a new
/// point to a linestring for example) and those retrieving a specific part
/// of a geometry (such as getting the 2nd interior ring of some polygon as a
/// linestring). However, this can easily be done by yourself by converting them
/// from/to coordinates with the `new_from_coordinates` and `to_coordinates`
/// methods.
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
    // TODO: reactivate when this binding is added in sfcgal-sys
    /*pub fn solid_set_exterior_shell(&self, shell: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Solid);
        precondition_match_type_other!(shell, GeomType::Polyhedralsurface);

        unsafe { sfcgal_solid_set_exterior_shell(self.c_geom.as_ptr(), shell.c_geom.as_ptr()) };

        Ok(())
    }*/

    /// Create a geometry by parsing a [WKT](https://en.wikipedia.org/wiki/Well-known_text) string.
    pub fn new(wkt: &str) -> Result<SFCGeometry> {
        initialize();

        let c_str = CString::new(wkt)?;

        let obj = unsafe { sfcgal_io_read_wkt(c_str.as_ptr(), wkt.len()) };

        unsafe { SFCGeometry::new_from_raw(obj, true) }
    }

    pub(crate) fn _init() {
        // Right now it is a no-op
        // Don't bother
        unsafe {
            sfcgal_init();
        }
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

    /// Returns a WKT representation of the given `SFCGeometry` using CGAL
    /// exact integer fractions as coordinate values.
    pub fn to_wkt(&self) -> Result<String> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_text(self.c_geom.as_ref(), ptr.as_mut_ptr(), &mut length);

            Ok(_c_string_with_size(ptr.assume_init(), length))
        }
    }

    /// Returns a WKT representation of the given `SFCGeometry` using
    /// floating point coordinate values with the desired number of
    /// decimals.
    pub fn to_wkt_decim(&self, nb_decim: i32) -> Result<String> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_text_decim(
                self.c_geom.as_ref(),
                nb_decim,
                ptr.as_mut_ptr(),
                &mut length,
            );

            Ok(_c_string_with_size(ptr.assume_init(), length))
        }
    }

    /// Creates a Wkb string of the given geometry. In memory version.
    pub fn to_wkb_in_memory(&self) -> Result<Vec<u8>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_wkb(self.c_geom.as_ref(), ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
    }

    /// # Safety
    /// Returns an EWKT representation of the given PreparedGeometry
    pub unsafe fn to_ewkt_in_memory(
        prepared: *const sfcgal_prepared_geometry_t,
        num_decimals: i32,
    ) -> Result<Vec<u8>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_prepared_geometry_as_ewkt(prepared, num_decimals, ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
    }

    /// Creates a Hexwkb string of the given geometry. In memory version.
    pub fn to_hexwkb_in_memory(&self) -> Result<Vec<u8>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_hexwkb(self.c_geom.as_ref(), ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
    }

    /// Read the binary prepared geometry
    pub fn io_read_binary_prepared(data: &[u8]) -> Result<*mut sfcgal_prepared_geometry_t> {
        unsafe {
            let ptr = data.as_ptr() as *const i8;

            let length = data.len();

            let result = sfcgal_io_read_binary_prepared(ptr, length);

            check_null_prepared_geom(result)?;

            Ok(result)
        }
    }

    /// Read the binary prepared geometry
    pub fn io_read_wkb(data: &[u8]) -> Result<*mut sfcgal_prepared_geometry_t> {
        unsafe {
            let ptr = data.as_ptr() as *const i8;

            let length = data.len();

            let result = sfcgal_io_read_wkb(ptr, length);

            check_null_prepared_geom(result)?;

            Ok(result)
        }
    }

    /// Read the binary prepared geometry
    pub fn io_read_ewkt(data: &[u8]) -> Result<*mut sfcgal_prepared_geometry_t> {
        unsafe {
            let ptr = data.as_ptr() as *const i8;

            let length = data.len();

            let result = sfcgal_io_read_ewkt(ptr, length);

            check_null_prepared_geom(result)?;

            Ok(result)
        }
    }

    /// # Safety
    /// Read the binary prepared geometry in memory
    pub unsafe fn io_write_binary_prepared(
        prepared_geometry: *mut sfcgal_prepared_geometry_t,
    ) -> Result<Vec<u8>> {
        unsafe {
            let mut ptr = MaybeUninit::<*mut c_char>::uninit();

            let mut length: usize = 0;

            sfcgal_io_write_binary_prepared(prepared_geometry, ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
    }

    /// Creates a OBJ file of the given geometry
    pub fn to_obj_file(&self, filename: &str) -> Result<()> {
        unsafe {
            let c_string = CString::new(filename)?;

            let raw: *mut c_char = c_string.into_raw();

            sfcgal_geometry_as_obj_file(self.c_geom.as_ptr(), raw);
        };

        Ok(())
    }

    /// Creates a OBJ string of the given geometry. In memory version.
    pub fn to_obj_in_memory(&self) -> Result<Vec<u8>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_obj(self.c_geom.as_ref(), ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
    }

    /// Creates a VTK file of the given geometry
    pub fn to_vtk_file(&self, filename: &str) -> Result<()> {
        unsafe {
            let c_string = CString::new(filename)?;

            let raw: *mut c_char = c_string.into_raw();

            sfcgal_geometry_as_vtk_file(self.c_geom.as_ptr(), raw);
        };

        Ok(())
    }

    /// Creates a VTK string of the given geometry. In memory version.
    pub fn to_vtk_in_memory(&self) -> Result<Vec<u8>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        let mut length: usize = 0;

        unsafe {
            sfcgal_geometry_as_vtk(self.c_geom.as_ref(), ptr.as_mut_ptr(), &mut length);

            Ok(get_raw_bytes(ptr.assume_init(), length))
        }
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

    /// Test if the given `SFCGeometry` is measured (has an 'm' coordinates)
    pub fn is_measured(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_measured(self.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Test if the given `SFCGeometry` is planar or not.
    pub fn is_planar(&self) -> Result<bool> {
        precondition_match_validity!(self);

        let rv = unsafe { sfcgal_geometry_is_planar(self.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Test if the given `SFCGeometry` is a 3d geometry or not.
    pub fn is_3d(&self) -> Result<bool> {
        let rv = unsafe { sfcgal_geometry_is_3d(self.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Returns reason for the invalidity or None in case of validity.
    pub fn validity_detail(&self) -> Result<Option<String>> {
        let mut ptr = MaybeUninit::<*mut c_char>::uninit();

        unsafe {
            let rv = sfcgal_geometry_is_valid_detail(
                self.c_geom.as_ptr(),
                ptr.as_mut_ptr(),
                std::ptr::null::<sfcgal_geometry_t>() as *mut *mut sfcgal_geometry_t,
            );

            match rv {
                1 => Ok(None),
                0 => Ok(Some(_string(ptr.assume_init()))),
                _ => Err(format_err!("SFCGAL error: {}", get_last_error())),
            }
        }
    }

    /// Returns the SFCGAL type of the given `SFCGeometry`.
    pub fn _type(&self) -> Result<GeomType> {
        let type_geom = unsafe { sfcgal_geometry_type_id(self.c_geom.as_ptr()) };

        GeomType::from_u32(type_geom)
            .ok_or_else(|| format_err!("Unknown geometry type (val={})", type_geom))
    }

    /// Computes the distance to an other `SFCGeometry`.
    pub fn distance(&self, other: &SFCGeometry) -> Result<f64> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let distance =
            unsafe { sfcgal_geometry_distance(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_computed_value(distance)
    }

    /// Computes the 3d distance to an other `SFCGeometry`.
    pub fn distance_3d(&self, other: &SFCGeometry) -> Result<f64> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let distance =
            unsafe { sfcgal_geometry_distance_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_computed_value(distance)
    }

    /// Computes the area of the given `SFCGeometry`.
    pub fn area(&self) -> Result<f64> {
        precondition_match_validity!(self);

        let area = unsafe { sfcgal_geometry_area(self.c_geom.as_ptr()) };

        check_computed_value(area)
    }

    /// Computes the 3d area of the given `SFCGeometry`.
    pub fn area_3d(&self) -> Result<f64> {
        precondition_match_validity!(self);

        let area = unsafe { sfcgal_geometry_area_3d(self.c_geom.as_ptr()) };

        check_computed_value(area)
    }

    /// Computes the volume of the given `SFCGeometry` (must be a volume).
    pub fn volume(&self) -> Result<f64> {
        precondition_match_validity!(self);

        let volume = unsafe { sfcgal_geometry_volume(self.c_geom.as_ptr()) };

        check_computed_value(volume)
    }

    /// Computes the orientation of the given `SFCGeometry` (must be a
    /// Polygon)
    pub fn orientation(&self) -> Result<Orientation> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_match_validity!(self);

        let orientation = unsafe { sfcgal_geometry_orientation(self.c_geom.as_ptr()) };

        Orientation::from_i32(orientation)
            .ok_or_else(|| format_err!("Error while retrieving orientation (val={})", orientation))
    }

    /// Test the intersection with an other `SFCGeometry`.
    pub fn intersects(&self, other: &SFCGeometry) -> Result<bool> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let rv = unsafe { sfcgal_geometry_intersects(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Test the 3d intersection with an other `SFCGeometry`.
    pub fn intersects_3d(&self, other: &SFCGeometry) -> Result<bool> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let rv =
            unsafe { sfcgal_geometry_intersects_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Returns the intersection of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn intersection(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_intersection(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d intersection of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn intersection_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_intersection_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Test if the geometry covers an other `SFCGeometry`.
    pub fn covers(&self, other: &SFCGeometry) -> Result<bool> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let rv = unsafe { sfcgal_geometry_covers(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Test if the geometry covers an other `SFCGeometry` in 3D.
    pub fn covers_3d(&self, other: &SFCGeometry) -> Result<bool> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let rv = unsafe { sfcgal_geometry_covers_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        check_predicate(rv)
    }

    /// Returns the difference of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn difference(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_difference(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d difference of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn difference_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_difference_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the union of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn union(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result = unsafe { sfcgal_geometry_union(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d union of the given `SFCGeometry` to an other
    /// `SFCGeometry`.
    pub fn union_3d(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_union_3d(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the minkowski sum of the given `SFCGeometry` and an other
    /// `SFCGEOMETRY`.
    pub fn minkowski_sum(&self, other: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        precondition_match_validity_other!(other);

        let result =
            unsafe { sfcgal_geometry_minkowski_sum(self.c_geom.as_ptr(), other.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the straight skeleton of the given `SFCGeometry`.
    pub fn straight_skeleton(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_straight_skeleton(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the straight skeleton of the given `SFCGeometry` with the
    /// distance to the border as M coordinate.
    pub fn straight_skeleton_distance_in_m(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result =
            unsafe { sfcgal_geometry_straight_skeleton_distance_in_m(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the extrude straight skeleton of the given Polygon.
    pub fn extrude_straight_skeleton(&self, height: f64) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_match_validity!(self);

        if height == 0. {
            bail!("Height cannot be 0.");
        }

        let result =
            unsafe { sfcgal_geometry_extrude_straight_skeleton(self.c_geom.as_ptr(), height) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the union of the polygon z-extrusion (with respect to
    /// building_height) and the extrude straight skeleton (with
    /// respect to roof_height) of the given Polygon.
    pub fn extrude_polygon_straight_skeleton(
        &self,
        building_height: f64,
        roof_height: f64,
    ) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_match_validity!(self);

        let result = unsafe {
            sfcgal_geometry_extrude_polygon_straight_skeleton(
                self.c_geom.as_ptr(),
                building_height,
                roof_height,
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the approximate medial axis for the given `SFCGeometry`
    /// Polygon.
    pub fn approximate_medial_axis(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_approximate_medial_axis(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the offset polygon of the given `SFCGeometry`.
    pub fn offset_polygon(&self, radius: f64) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_offset_polygon(self.c_geom.as_ptr(), radius) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the extrusion of the given `SFCGeometry` (not supported on
    /// Solid and Multisolid).
    pub fn extrude(&self, ex: f64, ey: f64, ez: f64) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_extrude(self.c_geom.as_ptr(), ex, ey, ez) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns a tesselation of the given `SFCGeometry`.
    pub fn tesselate(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_tesselate(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns a triangulation of the given `SFCGeometry`.
    pub fn triangulate_2dz(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_triangulate_2dz(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the convex hull of the given `SFCGeometry`.
    pub fn convexhull(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_convexhull(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the 3d convex hull of the given `SFCGeometry`.
    pub fn convexhull_3d(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_convexhull_3d(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the substring of the given `SFCGeometry` LineString between
    /// fractional distances.
    pub fn line_substring(&self, start: f64, end: f64) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Linestring);

        precondition_match_validity!(self);

        precondition_index_in_range!(start, (-1.0..=1.0));

        precondition_index_in_range!(end, (-1.0..=1.0));

        let result = unsafe { sfcgal_geometry_line_sub_string(self.c_geom.as_ptr(), start, end) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the alpha shape of the given `SFCGeometry` Point set.
    pub fn alpha_shapes(&self, alpha: f64, allow_holes: bool) -> Result<SFCGeometry> {
        if !self.is_valid().unwrap() {
            return Err(format_err!(
                "Error: alpha shapes can only be computed on valid geometries"
            ));
        }

        if alpha < 0.0 || !alpha.is_finite() {
            return Err(format_err!(
                "Error: alpha parameter must be positive or equal to 0.0, got {}",
                alpha,
            ));
        }

        let result =
            unsafe { sfcgal_geometry_alpha_shapes(self.c_geom.as_ptr(), alpha, allow_holes) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Return the optimal alpha shape of the given `SFCGeometry` Point set.
    pub fn optimal_alpha_shapes(
        &self,
        allow_holes: bool,
        nb_components: usize,
    ) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe {
            sfcgal_geometry_optimal_alpha_shapes(self.c_geom.as_ptr(), allow_holes, nb_components)
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Create a SFCGeometry collection type (MultiPoint, MultiLineString,
    /// MultiPolygon, MultiSolid or GeometryCollection) given a
    /// mutable slice of `SFCGeometry`'s (this is a destructive
    /// operation) ``` rust
    /// use sfcgal::SFCGeometry;
    /// let a = SFCGeometry::new("POINT (1.0 1.0)").unwrap();
    /// let b = SFCGeometry::new("POINT (2.0 2.0)").unwrap();
    /// let g = SFCGeometry::create_collection(&mut[a, b]).unwrap();
    /// assert_eq!(
    ///     g.to_wkt_decim(1).unwrap(),
    ///     "MULTIPOINT ((1.0 1.0),(2.0 2.0))",
    /// );
    /// ```
    pub fn create_collection(geoms: &mut [SFCGeometry]) -> Result<SFCGeometry> {
        if geoms.is_empty() {
            let res_geom = unsafe { sfcgal_geometry_collection_create() };

            return unsafe { SFCGeometry::new_from_raw(res_geom, true) };
        }

        let types = geoms
            .iter()
            .map(|g| g._type().unwrap())
            .collect::<Vec<GeomType>>();

        let multis = types
            .iter()
            .map(|gt| gt.is_collection_type())
            .collect::<Vec<bool>>();

        if !is_all_same(&types) || multis.iter().any(|&x| x) {
            let res_geom = unsafe { sfcgal_geometry_collection_create() };

            make_multi_geom(res_geom, geoms)
        } else if types[0] == GeomType::Point {
            let res_geom = unsafe { sfcgal_multi_point_create() };

            make_multi_geom(res_geom, geoms)
        } else if types[0] == GeomType::Linestring {
            let res_geom = unsafe { sfcgal_multi_linestring_create() };

            make_multi_geom(res_geom, geoms)
        } else if types[0] == GeomType::Polygon {
            let res_geom = unsafe { sfcgal_multi_polygon_create() };

            make_multi_geom(res_geom, geoms)
        } else if types[0] == GeomType::Solid {
            let mut res_geom = SFCGeometry::new("MULTISOLID EMPTY")?;

            res_geom.owned = false;

            make_multi_geom(res_geom.c_geom.as_ptr(), geoms)
        } else {
            unreachable!();
        }
    }

    /// Get the members of a SFCGeometry.
    /// Returns Err if the SFCGeometry if not a collection (i.e. if it's
    /// type is not in { MultiPoint, MultiLineString, MultiPolygon,
    /// MultiSolid, GeometryCollection }). The original geometry
    /// stay untouched. ``` rust
    /// use sfcgal::SFCGeometry;
    /// let g = SFCGeometry::new("MULTIPOINT ((1.0 1.0),(2.0
    /// 2.0))").unwrap(); let members =
    /// g.get_collection_members().unwrap(); assert_eq!(
    ///     members[0].to_wkt_decim(1).unwrap(),
    ///     "POINT (1.0 1.0)",
    /// );
    /// assert_eq!(
    ///     members[1].to_wkt_decim(1).unwrap(),
    ///     "POINT (2.0 2.0)",
    /// );
    /// ```
    pub fn get_collection_members(self) -> Result<Vec<SFCGeometry>> {
        let _type = self._type()?;

        if !_type.is_collection_type() {
            return Err(format_err!(
                "Error: the given geometry doesn't have any member ({:?} is not a collection type)",
                _type,
            ));
        }

        unsafe {
            let ptr = self.c_geom.as_ptr();

            let n_geom = sfcgal_geometry_collection_num_geometries(ptr);

            let mut result = Vec::new();

            for n in 0..n_geom {
                let _original_c_geom = sfcgal_geometry_collection_geometry_n(ptr, n);

                let clone_c_geom = sfcgal_geometry_clone(_original_c_geom);

                result.push(SFCGeometry::new_from_raw(clone_c_geom, true)?);
            }

            Ok(result)
        }
    }

    /// Creates an empty point
    pub fn point_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_point_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Creates a point from two X and Y coordinates
    pub fn point_create_from_xy(x: f64, y: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_point_create_from_xy(x, y) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Creates a point from three X, Y and M coordinates
    pub fn point_create_from_xym(x: f64, y: f64, m: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_point_create_from_xym(x, y, m) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Creates a point from three X, Y and Z coordinates
    pub fn point_create_from_xyz(x: f64, y: f64, z: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_point_create_from_xyz(x, y, z) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Create a point from x, y, z, m components.
    pub fn point_create_from_xyzm(&self, x: f64, y: f64, z: f64, m: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_point_create_from_xyzm(x, y, z, m) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the X coordinate of the given Point
    pub fn point_x(&self) -> Result<f64> {
        precondition_match_type!(self, GeomType::Point);

        precondition_match_not_empty!(self);

        unsafe { Ok(sfcgal_point_x(self.c_geom.as_ptr())) }
    }

    /// Returns the Y coordinate of the given Point
    pub fn point_y(&self) -> Result<f64> {
        precondition_match_type!(self, GeomType::Point);

        precondition_match_not_empty!(self);

        unsafe { Ok(sfcgal_point_y(self.c_geom.as_ptr())) }
    }

    /// Returns the Z coordinate of the given Point
    pub fn point_z(&self) -> Result<f64> {
        precondition_match_type!(self, GeomType::Point);

        precondition_match_not_empty!(self);

        unsafe {
            let result = sfcgal_point_z(self.c_geom.as_ptr());

            check_nan_value(result)
        }
    }

    /// Returns the M coordinate of the given Point
    pub fn point_m(&self) -> Result<f64> {
        precondition_match_type!(self, GeomType::Point);

        precondition_match_not_empty!(self);

        unsafe {
            let result = sfcgal_point_m(self.c_geom.as_ptr());

            check_nan_value(result)
        }
    }

    /// Sets one vertex of a Triangle
    pub fn triangle_set_vertex(&self, index: i32, vertex: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Triangle);

        precondition_index_in_range!(index, (0..3));

        unsafe {
            let explicit_converted_int: ::std::os::raw::c_int = index;

            sfcgal_triangle_set_vertex(
                self.c_geom.as_ptr(),
                explicit_converted_int,
                vertex.c_geom.as_ptr(),
            );

            Ok(())
        }
    }

    /// Sets one vertex of a Triangle from two coordinates
    pub fn triangle_set_vertex_from_xy(&self, index: i32, x: f64, y: f64) -> Result<()> {
        precondition_match_type!(self, GeomType::Triangle);

        precondition_index_in_range!(index, (0..3));

        unsafe {
            let explicit_converted_int: ::std::os::raw::c_int = index;

            sfcgal_triangle_set_vertex_from_xy(self.c_geom.as_ptr(), explicit_converted_int, x, y);

            Ok(())
        }
    }

    /// Sets one vertex of a Triangle from three coordinates
    pub fn triangle_set_vertex_from_xyz(&self, index: i32, x: f64, y: f64, z: f64) -> Result<()> {
        precondition_match_type!(self, GeomType::Triangle);

        precondition_index_in_range!(index, (0..3));

        unsafe {
            let explicit_converted_int: ::std::os::raw::c_int = index;

            sfcgal_triangle_set_vertex_from_xyz(
                self.c_geom.as_ptr(),
                explicit_converted_int,
                x,
                y,
                z,
            );

            Ok(())
        }
    }

    /// Returns the straight skeleton partition for the given Polygon
    pub fn straight_skeleton_partition(&self, auto_orientation: bool) -> Result<SFCGeometry> {
        match self._type()? {
            GeomType::Polygon | GeomType::Multipolygon | GeomType::Triangle => (),
            _ => {
                bail!("Wrong input Geometry type");
            }
        }

        precondition_match_validity!(self);

        let result = unsafe {
            sfcgal_geometry_straight_skeleton_partition(self.c_geom.as_ptr(), auto_orientation)
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the visibility polygon of a Point inside a Polygon
    pub fn visibility_point(&self, point: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_match_type_other!(point, GeomType::Point);

        precondition_match_validity!(self);

        if !self.covers(point)? {
            bail!("Point must be inside the Polygon");
        }

        let result = unsafe {
            sfcgal_geometry_visibility_point(self.c_geom.as_ptr(), point.c_geom.as_ptr())
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a geometry around the origin (0,0,0) by a given angle
    pub fn rotate(&self, angle: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_rotate(self.c_geom.as_ptr(), angle) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a geometry around the X axis by a given angle
    pub fn rotate_x(&self, angle: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_rotate_x(self.c_geom.as_ptr(), angle) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a geometry around the Y axis by a given angle
    pub fn rotate_y(&self, angle: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_rotate_y(self.c_geom.as_ptr(), angle) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a geometry around the Z axis by a given angle
    pub fn rotate_z(&self, angle: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_rotate_z(self.c_geom.as_ptr(), angle) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a geometry around a specified point by a given angle
    pub fn rotate_2d(&self, angle: f64, origin_x: f64, origin_y: f64) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_rotate_2d(self.c_geom.as_ptr(), angle, origin_x, origin_y) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a 3D geometry around a specified axis by a given angle
    pub fn rotate_3d(
        &self,
        angle: f64,
        axis_x_angle: f64,
        axis_y_angle: f64,
        axis_z_angle: f64,
    ) -> Result<SFCGeometry> {
        let result = unsafe {
            sfcgal_geometry_rotate_3d(
                self.c_geom.as_ptr(),
                angle,
                axis_x_angle,
                axis_y_angle,
                axis_z_angle,
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Rotates a 3D geometry around a specified axis and center point by a
    /// given
    #[allow(clippy::too_many_arguments)]
    pub fn rotate_3d_around_center(
        &self,
        angle: f64,
        axis_x_angle: f64,
        axis_y_angle: f64,
        axis_z_angle: f64,
        center_x: f64,
        center_y: f64,
        center_z: f64,
    ) -> Result<SFCGeometry> {
        let result = unsafe {
            sfcgal_geometry_rotate_3d_around_center(
                self.c_geom.as_ptr(),
                angle,
                axis_x_angle,
                axis_y_angle,
                axis_z_angle,
                center_x,
                center_y,
                center_z,
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Force a Right Handed Rule on the given Geometry
    pub fn force_rhr(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_force_rhr(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Force a Left Handed Rule on the given Geometry
    pub fn force_lhr(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_force_lhr(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Scale a geometry by a given factor
    pub fn scale(&self, scale: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_scale(self.c_geom.as_ptr(), scale) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Scale a geometry by different factors for each dimension
    pub fn scale_3d(&self, scale_x: f64, scale_y: f64, scale_z: f64) -> Result<SFCGeometry> {
        let result =
            unsafe { sfcgal_geometry_scale_3d(self.c_geom.as_ptr(), scale_x, scale_y, scale_z) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Scale a geometry by different factors for each dimension around a
    /// center
    pub fn scale_3d_around_center(
        &self,
        factor_x: f64,
        factor_y: f64,
        factor_z: f64,
        center_x: f64,
        center_y: f64,
        center_z: f64,
    ) -> Result<SFCGeometry> {
        let result = unsafe {
            sfcgal_geometry_scale_3d_around_center(
                self.c_geom.as_ptr(),
                factor_x,
                factor_y,
                factor_z,
                center_x,
                center_y,
                center_z,
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Round coordinates of the given Geometry
    pub fn round(&self, value: i32) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let explicit_conversion: ::std::os::raw::c_int = value;

        let result = unsafe { sfcgal_geometry_round(self.c_geom.as_ptr(), explicit_conversion) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Computes a 3D buffer around a geometry
    pub fn buffer3d(
        &self,
        radius: f64,
        segments: i32,
        buffer_type: BufferType,
    ) -> Result<SFCGeometry> {
        match self._type()? {
            GeomType::Point | GeomType::Linestring => (),
            _ => bail!("Geometry must be a Point or a Linestring"),
        }

        precondition_match_validity!(self);

        if radius <= 0. {
            bail!("Radius must be greater than 0.0");
        }

        if segments <= 3 {
            bail!("The algorithm needs at least 3 segments");
        }

        let explicit = {
            match buffer_type {
                BufferType::Round => 0_u32,
                BufferType::CylSphere => 1_u32,
                BufferType::Flat => 2_u32,
            }
        };

        unsafe {
            let result = sfcgal_geometry_buffer3d(self.c_geom.as_ptr(), radius, segments, explicit);

            SFCGeometry::new_from_raw(result, true)
        }
    }

    /// Gets the validity flag of the geometry
    pub fn has_validity_flag(&self) -> i32 {
        unsafe { sfcgal_geometry_has_validity_flag(self.c_geom.as_ptr()) }
    }

    /// Build the visibility polygon of the segment [pointA ; pointB] on a
    /// Polygon
    pub fn visibility_segment(
        &self,
        point_a: &SFCGeometry,
        point_b: &SFCGeometry,
    ) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe {
            sfcgal_geometry_visibility_segment(
                self.c_geom.as_ptr(),
                point_a.c_geom.as_ptr(),
                point_b.c_geom.as_ptr(),
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Convert a PolyhedralSurface to a Solid
    pub fn make_solid(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_geometry_make_solid(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Translate a geometry by a 2D vector
    pub fn translate_2d(&self, dx: f64, dy: f64) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_geometry_translate_2d(self.c_geom.as_ptr(), dx, dy) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Translate a geometry by a 3D vector
    pub fn translate_3d(
        &self,
        translation_x: f64,
        translation_y: f64,
        translation_z: f64,
    ) -> Result<SFCGeometry> {
        let result = unsafe {
            sfcgal_geometry_translate_3d(
                self.c_geom.as_ptr(),
                translation_x,
                translation_y,
                translation_z,
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Sets the validity flag of the geometry.
    pub fn force_valid(&self, validity: i32) {
        unsafe { sfcgal_geometry_force_valid(self.c_geom.as_ptr(), validity) };
    }

    /// Set the geometry validation mode
    pub fn set_geometry_validation(enabled: i32) {
        unsafe { sfcgal_set_geometry_validation(enabled) };
    }

    /// Creates an empty Triangle
    pub fn triangle_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_triangle_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns one the Triangle's vertex as a Point
    pub fn triangle_vertex(&self, index: i32) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Triangle);

        precondition_index_in_range!(index, (0..3));

        let result = unsafe { sfcgal_triangle_vertex(self.c_geom.as_ptr(), index) };

        let convert_mutability = result as *mut c_void;

        unsafe { SFCGeometry::new_from_raw(convert_mutability, true) }
    }

    /// Creates a Triangle from three given Point
    pub fn triangle_create_from_points(
        point_a: &SFCGeometry,
        point_b: &SFCGeometry,
        point_c: &SFCGeometry,
    ) -> Result<SFCGeometry> {
        precondition_match_type!(point_a, GeomType::Point);

        precondition_match_type!(point_a, GeomType::Point);

        precondition_match_type!(point_a, GeomType::Point);

        let result = unsafe {
            sfcgal_triangle_create_from_points(
                point_a.c_geom.as_ptr(),
                point_b.c_geom.as_ptr(),
                point_c.c_geom.as_ptr(),
            )
        };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Creates an empty TriangulatedSurface
    pub fn triangulated_surface_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_triangulated_surface_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the ith Triangle of a given TriangulatedSurface
    pub fn triangulated_surface_triangle_n(&self, index: usize) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Triangulatedsurface);

        precondition_index_in_result_value!(index, self.triangulated_surface_num_triangles());

        unsafe {
            let result = sfcgal_triangulated_surface_triangle_n(self.c_geom.as_ptr(), index);

            let convert_mutability = result as *mut c_void;

            SFCGeometry::new_from_raw(convert_mutability, true)
        }
    }

    /// Adds a Triangle to a given TriangulatedSurface
    pub fn triangulated_surface_add_triangle(&self, triangle: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Triangulatedsurface);

        precondition_match_type_other!(triangle, GeomType::Triangle);

        unsafe {
            sfcgal_triangulated_surface_add_triangle(self.c_geom.as_ptr(), triangle.c_geom.as_ptr())
        };

        Ok(())
    }

    /// Returns the number of triangles of a given TriangulatedSurface
    pub fn triangulated_surface_num_triangles(&self) -> Result<usize> {
        precondition_match_type!(self, GeomType::Triangulatedsurface);

        unsafe {
            Ok(sfcgal_triangulated_surface_num_triangles(
                self.c_geom.as_ptr(),
            ))
        }
    }

    /// Returns the optimal convex partition of a geometry (polygon without
    /// hole)
    pub fn optimal_convex_partition_2(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_optimal_convex_partition_2(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the approximal convex partition of a geometry (polygon
    /// without hole)
    pub fn approx_convex_partition_2(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_approx_convex_partition_2(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the greene approximal convex partition of a geometry
    /// (polygon without hole)
    pub fn greene_approx_convex_partition_2(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_greene_approx_convex_partition_2(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the y monotone partition of a geometry (polygon without
    /// hole)
    pub fn y_monotone_partition_2(&self) -> Result<SFCGeometry> {
        precondition_match_validity!(self);

        let result = unsafe { sfcgal_y_monotone_partition_2(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the number of points of the given LineString
    pub fn linestring_num_points(&self) -> Result<usize> {
        precondition_match_type!(self, GeomType::Linestring);

        unsafe { Ok(sfcgal_linestring_num_points(self.c_geom.as_ptr())) }
    }

    /// Adds a point to a LineString
    pub fn linestring_add_point(&self, point: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Linestring);

        precondition_match_type!(point, GeomType::Point);

        unsafe { sfcgal_linestring_add_point(self.c_geom.as_ptr(), point.c_geom.as_ptr()) };

        Ok(())
    }

    /// Returns the ith point of a given LineString
    pub fn linestring_point_n(&self, index: usize) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Point);

        match self.linestring_num_points() {
            Ok(num) => {
                if index >= num {
                    bail!("Overflowing index")
                }
            }
            Err(e) => bail!(e),
        }

        unsafe {
            let result = sfcgal_linestring_point_n(self.c_geom.as_ptr(), index);

            let convert_mutability = result as *mut c_void;

            SFCGeometry::new_from_raw(convert_mutability, true)
        }
    }

    /// Creates an empty LineString
    pub fn linestring_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_linestring_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Creates an empty PolyhedralSurface
    pub fn polyhedral_surface_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_polyhedral_surface_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the number of polygons of a given PolyhedralSurface
    pub fn polyhedral_surface_num_polygons(&self) -> Result<usize> {
        precondition_match_type!(self, GeomType::Polyhedralsurface);

        unsafe { Ok(sfcgal_polyhedral_surface_num_polygons(self.c_geom.as_ptr())) }
    }

    /// Adds a Polygon to a given PolyhedralSurface
    pub fn polyhedral_surface_add_polygon(&self, other: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Polyhedralsurface);

        precondition_match_type_other!(other, GeomType::Polygon);

        unsafe {
            sfcgal_polyhedral_surface_add_polygon(self.c_geom.as_ptr(), other.c_geom.as_ptr())
        };

        Ok(())
    }

    /// Returns the ith polygon of a given PolyhedralSurface
    pub fn polyhedral_surface_polygon_n(&self, index: usize) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Polyhedralsurface);

        precondition_index_in_result_value!(index, self.polyhedral_surface_num_polygons());

        let result = unsafe { sfcgal_polyhedral_surface_polygon_n(self.c_geom.as_ptr(), index) };

        let convert_mutability = result as *mut c_void;

        unsafe { SFCGeometry::new_from_raw(convert_mutability, true) }
    }

    /// Sets the error handlers. These callbacks are called on warning or
    /// error
    pub fn set_error_handlers(
        warning_handler: sfcgal_error_handler_t,
        error_handler: sfcgal_error_handler_t,
    ) {
        unsafe { sfcgal_set_error_handlers(warning_handler, error_handler) }
    }

    /// Sets the error handlers. These callbacks are called on warning or
    /// error
    pub fn set_alloc_handlers(
        alloc_handler: sfcgal_alloc_handler_t,
        free_handler: sfcgal_free_handler_t,
    ) {
        unsafe { sfcgal_set_alloc_handlers(alloc_handler, free_handler) };
    }

    /// Creates an empty Solid
    pub fn solid_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_solid_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the number of shells of a given Solid
    pub fn solid_num_shells(&self) -> Result<usize> {
        precondition_match_type!(self, GeomType::Solid);

        unsafe { Ok(sfcgal_solid_num_shells(self.c_geom.as_ptr())) }
    }

    /// Creates a Solid from an exterior shell
    pub fn solid_create_from_exterior_shell(exterior_shell: &SFCGeometry) -> Result<SFCGeometry> {
        precondition_match_type_other!(exterior_shell, GeomType::Polyhedralsurface);

        unsafe {
            let result = sfcgal_solid_create_from_exterior_shell(exterior_shell.c_geom.as_ptr());

            SFCGeometry::new_from_raw(result, true)
        }
    }

    /// Returns the ith shell of a given Solid
    pub fn solid_shell_n(&self, index: usize) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Solid);

        precondition_index_in_result_value!(index, self.solid_num_shells());

        let result = unsafe { sfcgal_solid_shell_n(self.c_geom.as_ptr(), index) };

        let convert_mutability = result as *mut c_void;

        unsafe { SFCGeometry::new_from_raw(convert_mutability, true) }
    }

    /// Adds a shell to a given Solid
    pub fn solid_add_interior_shell(&self, shell: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Solid);

        precondition_match_type_other!(shell, GeomType::Polyhedralsurface);

        unsafe { sfcgal_solid_add_interior_shell(self.c_geom.as_ptr(), shell.c_geom.as_ptr()) };

        Ok(())
    }

    /// Creates an empty Polygon
    pub fn polygon_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_polygon_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Returns the number of interior rings of a given Polygon
    pub fn polygon_num_interior_rings(&self) -> Result<usize> {
        precondition_match_type!(self, GeomType::Polygon);

        unsafe { Ok(sfcgal_polygon_num_interior_rings(self.c_geom.as_ptr())) }
    }

    /// Returns the exterior ring of a given Polygon
    pub fn polygon_exterior_ring(&self) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_polygon_exterior_ring(self.c_geom.as_ptr()) };

        let convert_mutability = result as *mut c_void;

        unsafe { SFCGeometry::new_from_raw(convert_mutability, true) }
    }

    /// Creates an empty Polygon from an extrior ring
    pub fn polygon_create_from_exterior_ring(&self) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Linestring);

        let result = unsafe { sfcgal_polygon_create_from_exterior_ring(self.c_geom.as_ptr()) };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    /// Adds an interior ring to a given Polygon
    pub fn polygon_add_interior_ring(&self, ring: &SFCGeometry) -> Result<()> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_match_type!(ring, GeomType::Linestring);

        unsafe { sfcgal_polygon_add_interior_ring(self.c_geom.as_ptr(), ring.c_geom.as_ptr()) };

        Ok(())
    }

    /// Returns the ith interior ring of a given Polygon
    pub fn polygon_interior_ring_n(&self, index: usize) -> Result<SFCGeometry> {
        precondition_match_type!(self, GeomType::Polygon);

        precondition_index_in_result_value!(index, self.polygon_num_interior_rings());

        let result = unsafe { sfcgal_polygon_interior_ring_n(self.c_geom.as_ptr(), index) };

        let convert_mutability = result as *mut c_void;

        unsafe { SFCGeometry::new_from_raw(convert_mutability, true) }
    }

    /// Creates an empty MultiSolid
    pub fn multi_solid_create() -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_multi_solid_create() };

        unsafe { SFCGeometry::new_from_raw(result, true) }
    }

    // ------------------------------
    // Prepared geometries section
    // ------------------------------

    /// # Safety
    /// Sets SRID associated with a given PreparedGeometry
    pub unsafe fn prepared_geometry_set_srid(
        prepared_geometry: *mut sfcgal_prepared_geometry_t,
        srid: srid_t,
    ) {
        sfcgal_prepared_geometry_set_srid(prepared_geometry, srid);
    }

    /// # Safety
    /// Deletes a given PreparedGeometry
    pub unsafe fn prepared_geometry_delete(prepared_geometry: *mut sfcgal_prepared_geometry_t) {
        sfcgal_prepared_geometry_delete(prepared_geometry);
    }

    /// # Safety
    /// Creates a PreparedGeometry from a Geometry and an SRID
    pub fn prepared_geometry_create_from_geometry(
        &self,
        srid: srid_t,
    ) -> Result<*mut sfcgal_prepared_geometry_t> {
        let result =
            unsafe { sfcgal_prepared_geometry_create_from_geometry(self.c_geom.as_ptr(), srid) };

        check_null_prepared_geom(result)?;

        Ok(result)
    }

    /// # Safety
    /// Sets the Geometry associated with the given PreparedGeometry
    pub unsafe fn prepared_geometry_set_geometry(&self, prepared: *mut sfcgal_prepared_geometry_t) {
        sfcgal_prepared_geometry_set_geometry(prepared, self.c_geom.as_ptr());
    }

    /// # Safety
    /// Returns SRID associated with a given PreparedGeometry
    pub unsafe fn prepared_geometry_srid(prepared: *mut sfcgal_prepared_geometry_t) -> u32 {
        sfcgal_prepared_geometry_srid(prepared)
    }

    /// Creates an empty PreparedGeometry
    pub fn prepared_geometry_create() -> Result<*mut sfcgal_prepared_geometry_t> {
        let result = unsafe { sfcgal_prepared_geometry_create() };

        check_null_prepared_geom(result)?;

        Ok(result)
    }

    /// # Safety
    /// Returns the Geometry associated with a given PreparedGeometry
    pub unsafe fn prepared_geometry_geometry(
        prepared_geometry: *mut sfcgal_prepared_geometry_t,
    ) -> Result<SFCGeometry> {
        let result = unsafe { sfcgal_prepared_geometry_geometry(prepared_geometry) };

        // NOTE: Not proud..
        let converted = result as *mut c_void;

        SFCGeometry::new_from_raw(converted, true)
    }
}
fn is_all_same<T>(arr: &[T]) -> bool
where
    T: Ord + Eq,
{
    arr.iter().min() == arr.iter().max()
}
fn make_multi_geom(
    out_multi: *mut sfcgal_geometry_t,
    geoms: &mut [SFCGeometry],
) -> Result<SFCGeometry> {
    for sfcgal_geom in geoms.iter_mut() {
        unsafe {
            sfcgal_geom.owned = false;

            sfcgal_geometry_collection_add_geometry(
                out_multi,
                sfcgal_geom.c_geom.as_ptr() as *mut sfcgal_geometry_t,
            )
        };
    }

    unsafe { SFCGeometry::new_from_raw(out_multi, true) }
}

pub fn _sfcgal_get_full_version() -> String {
    let result = unsafe { sfcgal_full_version() };

    _string(result)
}

pub fn _sfcgal_get_version() -> String {
    let result = unsafe { sfcgal_version() };

    _string(result)
}

#[cfg(test)]

mod tests {

    use std::{env, f64::consts::PI};

    use num_traits::abs;

    use super::*;
    use crate::{Point2d, Point3d, ToCoordinates};

    #[test]

    fn creation_point_from_wkt() {
        let geom = SFCGeometry::new("POINT (1.0 1.0)");

        assert!(geom.is_ok());
    }

    #[test]

    fn creation_polygon_from_wkt() {
        let geom = SFCGeometry::new("POLYGON((0.0 0.0, 1.0 0.0, 1.0 1.0, 0.0 0.0))");

        assert!(geom.is_ok());

        let geom = geom.unwrap();

        assert!(geom.is_valid().unwrap());

        let geom1 = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        assert!(geom.intersects(&geom1).unwrap());
    }

    #[test]

    fn writing_to_wkt() {
        let geom = SFCGeometry::new("POINT (1.0 1.0)");

        assert!(geom.is_ok());

        let wkt = geom.unwrap().to_wkt();

        assert!(wkt.is_ok());

        assert_eq!(wkt.unwrap(), String::from("POINT (1/1 1/1)"));
    }

    #[test]

    fn writing_to_wkt_with_decimals() {
        let geom = SFCGeometry::new("POINT (1.0 1.0)");

        assert!(geom.is_ok());

        let wkt = geom.unwrap().to_wkt_decim(1);

        assert!(wkt.is_ok());

        assert_eq!(wkt.unwrap(), String::from("POINT (1.0 1.0)"));
    }

    #[test]

    fn creation_failed_with_error_message() {
        let geom = SFCGeometry::new("POINT (1, 1)");

        assert!(geom.is_err());

        assert_eq!(geom.err().unwrap().to_string(), "Obtained null pointer when creating geometry: WKT parse error, Coordinate dimension < 2 (, 1))",)
    }

    #[test]

    fn distance_to_other() {
        let pt1 = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        let pt2 = SFCGeometry::new("POINT (10.0 1.0)").unwrap();

        let distance = pt1.distance(&pt2).unwrap();

        assert_eq!(distance, 9.0);
    }

    #[test]

    fn distance_3d_to_other() {
        let pt1 = SFCGeometry::new("POINT (1.0 1.0 2.0)").unwrap();

        let pt2 = SFCGeometry::new("POINT (10.0 1.0 2.0)").unwrap();

        let distance = pt1.distance_3d(&pt2).unwrap();

        assert_eq!(distance, 9.0);
    }

    #[test]

    fn measured_geometry() {
        let pt1 = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        let pt2 = SFCGeometry::new("POINTM(1.0 1.0 2.0)").unwrap();

        assert!(!pt1.is_measured().unwrap());

        assert!(pt2.is_measured().unwrap());
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

        assert!(surface.volume().is_err());
    }

    #[test]

    fn predicates() {
        let pt = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        assert!(pt.is_valid().unwrap());

        assert!(!pt.is_3d().unwrap());

        assert!(!pt.is_empty().unwrap());

        assert_eq!(
            pt.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let linestring_3d = SFCGeometry::new("LINESTRING (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        assert!(linestring_3d.is_valid().unwrap());

        assert!(linestring_3d.is_3d().unwrap());

        assert!(!linestring_3d.is_empty().unwrap());

        assert_eq!(
            linestring_3d.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let empty_geom = SFCGeometry::new("LINESTRING EMPTY").unwrap();

        assert!(empty_geom.is_valid().unwrap());

        assert!(!empty_geom.is_3d().unwrap());

        assert!(empty_geom.is_empty().unwrap());

        assert_eq!(
            linestring_3d.is_planar().err().unwrap().to_string(),
            "SFCGAL error: is_planar() only applies to polygons",
        );

        let polyg = SFCGeometry::new("POLYGON ((1 1, 3 1, 4 4, 1 3, 1 1))").unwrap();

        assert!(polyg.is_valid().unwrap());

        assert!(!polyg.is_3d().unwrap());

        assert!(!polyg.is_empty().unwrap());

        assert!(polyg.is_planar().unwrap());

        assert!(pt.intersects(&polyg).unwrap());

        assert!(!pt.intersects_3d(&linestring_3d).unwrap());
    }

    #[test]

    fn validity_detail_on_valid_geom() {
        let line = SFCGeometry::new("LINESTRING (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        assert!(line.is_valid().unwrap());

        assert_eq!(line.validity_detail().unwrap(), None);
    }

    #[test]

    fn validity_detail_on_invalid_geom_1() {
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

        assert!(!surface.is_valid().unwrap());

        assert_eq!(surface.validity_detail().unwrap(), Some(String::from("inconsistent orientation of PolyhedralSurface detected at edge 3 (4-7) of polygon 5")),);
    }

    #[test]

    fn validity_detail_on_invalid_geom_2() {
        let surface = SFCGeometry::new("POLYGON ((1 2,1 2,1 2,1 2))").unwrap();

        assert!(!surface.is_valid().unwrap());

        assert_eq!(
            surface.validity_detail().unwrap(),
            Some(String::from("ring 0 degenerated to a point")),
        );
    }

    #[test]

    fn validity_detail_on_invalid_geom_3() {
        let surface = SFCGeometry::new("LINESTRING (1 2, 1 2, 1 2)").unwrap();

        assert!(!surface.is_valid().unwrap());

        assert_eq!(
            surface.validity_detail().unwrap(),
            Some(String::from("no length")),
        );
    }

    #[test]

    fn straight_skeleton() {
        let geom = SFCGeometry::new("POLYGON ((0 0,1 0,1 1,0 1,0 0))").unwrap();

        let result = geom.straight_skeleton().unwrap();

        let wkt = result.to_wkt_decim(1).unwrap();

        assert_eq!(wkt, "MULTILINESTRING ((0.0 0.0,0.5 0.5),(1.0 0.0,0.5 0.5),(1.0 1.0,0.5 0.5),(0.0 1.0,0.5 0.5))",);
    }

    #[test]

    fn straight_skeleton_distance_in_m() {
        let geom = SFCGeometry::new("POLYGON ((0 0,1 0,1 1,0 1,0 0))").unwrap();

        let result = geom.straight_skeleton_distance_in_m().unwrap();

        let wkt = result.to_wkt_decim(1).unwrap();

        assert_eq!(
            wkt,
            "MULTILINESTRING M (\
             (0.0 0.0 0.0,0.5 0.5 0.5),\
             (1.0 0.0 0.0,0.5 0.5 0.5),\
             (1.0 1.0 0.0,0.5 0.5 0.5),\
             (0.0 1.0 0.0,0.5 0.5 0.5))",
        );
    }

    #[test]

    fn extrude_straight_skeleton() {
        // Test adapted from https://gitlab.com/sfcgal/SFCGAL/-/blob/master/test/unit/SFCGAL/algorithm/StraightSkeletonTest.cpp#L275
        let geom = SFCGeometry::new("POLYGON ((0 0, 5 0, 5 5, 4 5, 4 4, 0 4, 0 0))").unwrap();

        let result = geom.extrude_straight_skeleton(2.).unwrap();

        let wkt = result.to_wkt_decim(2).unwrap();

        assert_eq!(
            wkt,
            "POLYHEDRALSURFACE Z (((4.00 5.00 0.00,5.00 5.00 0.00,4.00 4.00 0.00,4.00 \
              5.00 0.00)),((0.00 4.00 0.00,4.00 4.00 0.00,0.00 0.00 0.00,0.00 4.00 \
              0.00)),((4.00 4.00 0.00,5.00 0.00 0.00,0.00 0.00 0.00,4.00 4.00 \
              0.00)),((5.00 5.00 0.00,5.00 0.00 0.00,4.00 4.00 0.00,5.00 5.00 \
              0.00)),((0.00 4.00 0.00,0.00 0.00 0.00,2.00 2.00 2.00,0.00 4.00 \
              0.00)),((0.00 0.00 0.00,5.00 0.00 0.00,3.00 2.00 2.00,0.00 0.00 \
              0.00)),((2.00 2.00 2.00,0.00 0.00 0.00,3.00 2.00 2.00,2.00 2.00 \
              2.00)),((4.50 3.50 0.50,5.00 5.00 0.00,4.50 4.50 0.50,4.50 3.50 \
              0.50)),((3.00 2.00 2.00,5.00 0.00 0.00,4.50 3.50 0.50,3.00 2.00 \
              2.00)),((4.50 3.50 0.50,5.00 0.00 0.00,5.00 5.00 0.00,4.50 3.50 \
              0.50)),((5.00 5.00 0.00,4.00 5.00 0.00,4.50 4.50 0.50,5.00 5.00 \
              0.00)),((4.50 4.50 0.50,4.00 4.00 0.00,4.50 3.50 0.50,4.50 4.50 \
              0.50)),((4.50 4.50 0.50,4.00 5.00 0.00,4.00 4.00 0.00,4.50 4.50 \
              0.50)),((4.00 4.00 0.00,0.00 4.00 0.00,2.00 2.00 2.00,4.00 4.00 \
              0.00)),((4.50 3.50 0.50,4.00 4.00 0.00,3.00 2.00 2.00,4.50 3.50 \
              0.50)),((3.00 2.00 2.00,4.00 4.00 0.00,2.00 2.00 2.00,3.00 2.00 \
              2.00)))"
        );
    }

    #[test]

    fn extrude_polygon_straight_skeleton() {
        // Test adapted from https://gitlab.com/sfcgal/SFCGAL/-/blob/master/test/unit/SFCGAL/algorithm/StraightSkeletonTest.cpp#L354
        let geom = SFCGeometry::new(
            "POLYGON (( 0 0, 5 0, 5 5, 4 5, 4 4, 0 4, 0 0 ), (1 1, 1 2, 2 2, 2 1, 1 1))",
        )
        .unwrap();

        let result = geom.extrude_polygon_straight_skeleton(9., 2.).unwrap();

        let wkt = result.to_wkt_decim(1).unwrap();

        assert_eq!(
            wkt,
            "POLYHEDRALSURFACE Z (((0.0 0.0 0.0,0.0 4.0 0.0,4.0 4.0 \
             0.0,4.0 5.0 0.0,5.0 5.0 0.0,5.0 0.0 0.0,0.0 0.0 0.0),\
             (1.0 1.0 0.0,2.0 1.0 0.0,2.0 2.0 0.0,1.0 2.0 0.0,1.0 1.0 0.0)),\
             ((0.0 0.0 0.0,0.0 0.0 9.0,0.0 4.0 9.0,0.0 4.0 0.0,0.0 0.0 0.0)),\
             ((0.0 4.0 0.0,0.0 4.0 9.0,4.0 4.0 9.0,4.0 4.0 0.0,0.0 4.0 0.0)),\
             ((4.0 4.0 0.0,4.0 4.0 9.0,4.0 5.0 9.0,4.0 5.0 0.0,4.0 4.0 0.0)),\
             ((4.0 5.0 0.0,4.0 5.0 9.0,5.0 5.0 9.0,5.0 5.0 0.0,4.0 5.0 0.0)),\
             ((5.0 5.0 0.0,5.0 5.0 9.0,5.0 0.0 9.0,5.0 0.0 0.0,5.0 5.0 0.0)),\
             ((5.0 0.0 0.0,5.0 0.0 9.0,0.0 0.0 9.0,0.0 0.0 0.0,5.0 0.0 0.0)),\
             ((1.0 1.0 0.0,1.0 1.0 9.0,2.0 1.0 9.0,2.0 1.0 0.0,1.0 1.0 0.0)),\
             ((2.0 1.0 0.0,2.0 1.0 9.0,2.0 2.0 9.0,2.0 2.0 0.0,2.0 1.0 0.0)),\
             ((2.0 2.0 0.0,2.0 2.0 9.0,1.0 2.0 9.0,1.0 2.0 0.0,2.0 2.0 0.0)),\
             ((1.0 2.0 0.0,1.0 2.0 9.0,1.0 1.0 9.0,1.0 1.0 0.0,1.0 2.0 0.0)),\
             ((4.0 5.0 9.0,5.0 5.0 9.0,4.0 4.0 9.0,4.0 5.0 9.0)),\
             ((2.0 1.0 9.0,5.0 0.0 9.0,0.0 0.0 9.0,2.0 1.0 9.0)),\
             ((5.0 5.0 9.0,5.0 0.0 9.0,4.0 4.0 9.0,5.0 5.0 9.0)),\
             ((2.0 1.0 9.0,0.0 0.0 9.0,1.0 1.0 9.0,2.0 1.0 9.0)),\
             ((1.0 2.0 9.0,1.0 1.0 9.0,0.0 0.0 9.0,1.0 2.0 9.0)),\
             ((0.0 4.0 9.0,2.0 2.0 9.0,1.0 2.0 9.0,0.0 4.0 9.0)),\
             ((0.0 4.0 9.0,1.0 2.0 9.0,0.0 0.0 9.0,0.0 4.0 9.0)),\
             ((4.0 4.0 9.0,5.0 0.0 9.0,2.0 2.0 9.0,4.0 4.0 9.0)),\
             ((4.0 4.0 9.0,2.0 2.0 9.0,0.0 4.0 9.0,4.0 4.0 9.0)),\
             ((2.0 2.0 9.0,5.0 0.0 9.0,2.0 1.0 9.0,2.0 2.0 9.0)),\
             ((0.5 2.5 9.5,0.0 0.0 9.0,0.5 0.5 9.5,0.5 2.5 9.5)),\
             ((1.0 3.0 10.0,0.0 4.0 9.0,0.5 2.5 9.5,1.0 3.0 10.0)),\
             ((0.5 2.5 9.5,0.0 4.0 9.0,0.0 0.0 9.0,0.5 2.5 9.5)),\
             ((2.5 0.5 9.5,5.0 0.0 9.0,3.5 1.5 10.5,2.5 0.5 9.5)),\
             ((0.0 0.0 9.0,5.0 0.0 9.0,2.5 0.5 9.5,0.0 0.0 9.0)),\
             ((0.5 0.5 9.5,0.0 0.0 9.0,2.5 0.5 9.5,0.5 0.5 9.5)),\
             ((4.5 3.5 9.5,5.0 5.0 9.0,4.5 4.5 9.5,4.5 3.5 9.5)),\
             ((3.5 2.5 10.5,3.5 1.5 10.5,4.5 3.5 9.5,3.5 2.5 10.5)),\
             ((4.5 3.5 9.5,5.0 0.0 9.0,5.0 5.0 9.0,4.5 3.5 9.5)),\
             ((3.5 1.5 10.5,5.0 0.0 9.0,4.5 3.5 9.5,3.5 1.5 10.5)),\
             ((5.0 5.0 9.0,4.0 5.0 9.0,4.5 4.5 9.5,5.0 5.0 9.0)),\
             ((4.5 4.5 9.5,4.0 4.0 9.0,4.5 3.5 9.5,4.5 4.5 9.5)),\
             ((4.5 4.5 9.5,4.0 5.0 9.0,4.0 4.0 9.0,4.5 4.5 9.5)),\
             ((3.0 3.0 10.0,0.0 4.0 9.0,1.0 3.0 10.0,3.0 3.0 10.0)),\
             ((3.5 2.5 10.5,4.5 3.5 9.5,3.0 3.0 10.0,3.5 2.5 10.5)),\
             ((3.0 3.0 10.0,4.0 4.0 9.0,0.0 4.0 9.0,3.0 3.0 10.0)),\
             ((4.5 3.5 9.5,4.0 4.0 9.0,3.0 3.0 10.0,4.5 3.5 9.5)),\
             ((2.0 1.0 9.0,1.0 1.0 9.0,0.5 0.5 9.5,2.0 1.0 9.0)),\
             ((2.5 0.5 9.5,2.0 1.0 9.0,0.5 0.5 9.5,2.5 0.5 9.5)),\
             ((1.0 1.0 9.0,1.0 2.0 9.0,0.5 2.5 9.5,1.0 1.0 9.0)),\
             ((0.5 0.5 9.5,1.0 1.0 9.0,0.5 2.5 9.5,0.5 0.5 9.5)),\
             ((1.0 3.0 10.0,2.0 2.0 9.0,3.0 3.0 10.0,1.0 3.0 10.0)),\
             ((0.5 2.5 9.5,1.0 2.0 9.0,1.0 3.0 10.0,0.5 2.5 9.5)),\
             ((1.0 3.0 10.0,1.0 2.0 9.0,2.0 2.0 9.0,1.0 3.0 10.0)),\
             ((2.0 2.0 9.0,2.0 1.0 9.0,2.5 0.5 9.5,2.0 2.0 9.0)),\
             ((3.5 2.5 10.5,3.0 3.0 10.0,3.5 1.5 10.5,3.5 2.5 10.5)),\
             ((3.5 1.5 10.5,2.0 2.0 9.0,2.5 0.5 9.5,3.5 1.5 10.5)),\
             ((3.0 3.0 10.0,2.0 2.0 9.0,3.5 1.5 10.5,3.0 3.0 10.0)))"
        );
    }

    #[test]

    fn tesselate() {
        let geom = SFCGeometry::new("POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();

        let result = geom.tesselate().unwrap();

        let output_wkt = result.to_wkt_decim(1).unwrap();

        assert_eq!(
            output_wkt,
            "TIN (((0.0 1.0,1.0 0.0,1.0 1.0,0.0 1.0)),((0.0 1.0,0.0 0.0,1.0 0.0,0.0 1.0)))",
        );
    }

    #[test]

    fn offset_polygon() {
        let geom = SFCGeometry::new("POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();

        let buff = geom.offset_polygon(1.).unwrap();

        assert!(buff.is_valid().unwrap());

        assert!(!buff.is_empty().unwrap());
    }

    #[test]

    fn extrude_polygon() {
        let geom = SFCGeometry::new("POLYGON ((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0))").unwrap();

        let extr = geom.extrude(0., 0., 1.).unwrap();

        assert!(extr.is_valid().unwrap());

        assert!(!extr.is_empty().unwrap());

        assert_eq!(extr._type().unwrap(), GeomType::Solid);
    }

    #[test]

    fn tesselate_invariant_geom() {
        let input_wkt = String::from("POINT (1.0 1.0)");

        let pt = SFCGeometry::new(&input_wkt).unwrap();

        let result = pt.tesselate().unwrap();

        let output_wkt = result.to_wkt_decim(1).unwrap();

        assert_eq!(input_wkt, output_wkt);
    }

    #[test]

    fn line_substring() {
        let g = SFCGeometry::new("LINESTRING Z (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        let result = g.line_substring(-0.2, 0.2).unwrap();

        assert_eq!(
            result.to_wkt_decim(1).unwrap(),
            "LINESTRING Z (2.8 1.8 1.8,8.2 1.2 1.9)"
        );

        // With "start" or "end" point not in [-1; 1]
        assert_eq!(
            g.line_substring(-2., 0.2).err().unwrap().to_string(),
            "Index not in the expected range"
        );
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

        assert!(diff.is_valid().unwrap());

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

        assert!(diff.is_valid().unwrap());

        assert_ulps_eq!(diff.volume().unwrap(), 0.5);
    }

    #[test]

    fn alpha_shapes_on_point() {
        let multipoint = SFCGeometry::new("POINT (1 3)").unwrap();

        let res = multipoint.alpha_shapes(20.0, false).unwrap();

        assert_eq!(res.to_wkt_decim(1).unwrap(), "GEOMETRYCOLLECTION EMPTY");
    }

    #[test]

    fn alpha_shapes_on_multipoint() {
        let multipoint = SFCGeometry::new("MULTIPOINT ((1 2),(2 2),(3 0),(1 3))").unwrap();

        let res = multipoint.alpha_shapes(20.0, false).unwrap();

        assert_eq!(
            res.to_wkt_decim(1).unwrap(),
            "POLYGON ((1.0 2.0,1.0 3.0,2.0 2.0,3.0 0.0,1.0 2.0))"
        );
    }

    #[test]

    fn alpha_shapes_on_linestring() {
        let multipoint = SFCGeometry::new("LINESTRING (1 2,2 2,3 0,1 3)").unwrap();

        let res = multipoint.alpha_shapes(20.0, false).unwrap();

        assert_eq!(
            res.to_wkt_decim(1).unwrap(),
            "POLYGON ((1.0 2.0,1.0 3.0,2.0 2.0,3.0 0.0,1.0 2.0))"
        );
    }

    #[test]

    fn alpha_shapes_on_multilinestring() {
        let multipoint =
            SFCGeometry::new("MULTILINESTRING ((1 2,2 2,3 0,1 3), (2 6,3 5,4 2))").unwrap();

        let res = multipoint.alpha_shapes(20.0, false).unwrap();

        assert_eq!(
            res.to_wkt_decim(1).unwrap(),
            "POLYGON ((1.0 2.0,1.0 3.0,2.0 6.0,3.0 5.0,4.0 2.0,3.0 0.0,1.0 2.0))"
        );
    }

    #[test]

    fn alpha_shapes_on_polygon() {
        let pol1 = SFCGeometry::new("POLYGON ((0 0, 0 4, 4 4, 4 0, 0 0))").unwrap();

        let res = pol1.alpha_shapes(10.0, false).unwrap();

        assert_eq!(
            res.to_wkt_decim(1).unwrap(),
            "POLYGON ((0.0 0.0,0.0 4.0,4.0 4.0,4.0 0.0,0.0 0.0))"
        );
    }

    #[test]

    fn alpha_shapes_on_invalid() {
        let pol1 = SFCGeometry::new("LINESTRING (1 2, 1 2, 1 2, 1 2)").unwrap();

        let res = pol1.alpha_shapes(10.0, false);

        assert!(res.is_err());
    }

    #[test]

    fn create_collection_empty() {
        let g = SFCGeometry::create_collection(&mut []).unwrap();

        assert_eq!(g.to_wkt_decim(1).unwrap(), "GEOMETRYCOLLECTION EMPTY",);
    }

    #[test]

    fn create_collection_heterogenous() {
        let a = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        let b = SFCGeometry::new("LINESTRING Z (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        let g = SFCGeometry::create_collection(&mut [a, b]).unwrap();

        assert_eq!(
            g.to_wkt_decim(1).unwrap(),
            "GEOMETRYCOLLECTION (POINT (1.0 1.0),LINESTRING Z (10.0 1.0 2.0,1.0 2.0 1.7))",
        );
    }

    #[test]

    fn create_collection_multipoint_from_points() {
        let a = SFCGeometry::new("POINT (1.0 1.0)").unwrap();

        let b = SFCGeometry::new("POINT (2.0 2.0)").unwrap();

        let g = SFCGeometry::create_collection(&mut [a, b]).unwrap();

        assert_eq!(
            g.to_wkt_decim(1).unwrap(),
            "MULTIPOINT ((1.0 1.0),(2.0 2.0))",
        );
    }

    #[test]

    fn create_collection_multilinestring_from_linestrings() {
        let a = SFCGeometry::new("LINESTRING (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        let b = SFCGeometry::new("LINESTRING (10.0 1.0 2.0, 1.0 2.0 1.7)").unwrap();

        let g = SFCGeometry::create_collection(&mut [a, b]).unwrap();

        assert_eq!(
            g.to_wkt_decim(1).unwrap(),
            "MULTILINESTRING Z ((10.0 1.0 2.0,1.0 2.0 1.7),(10.0 1.0 2.0,1.0 2.0 1.7))",
        );
    }

    #[test]

    fn create_collection_multisolid_from_solids() {
        let a = SFCGeometry::new(
            "SOLID((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
             ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
             ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
             ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
             ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
             ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))",
        )
        .unwrap();

        let b = SFCGeometry::new(
            "SOLID((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
             ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
             ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
             ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
             ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
             ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))",
        )
        .unwrap();

        let g = SFCGeometry::create_collection(&mut [a, b]).unwrap();

        assert_eq!(
            g.to_wkt_decim(1).unwrap(),
            "MULTISOLID Z (\
             ((\
             ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 0.0,0.0 0.0 0.0)),\
             ((0.0 0.0 0.0,0.0 1.0 0.0,1.0 1.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0)),\
             ((0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0)),\
             ((1.0 0.0 0.0,1.0 1.0 0.0,1.0 1.0 1.0,1.0 0.0 1.0,1.0 0.0 0.0)),\
             ((0.0 0.0 1.0,1.0 0.0 1.0,1.0 1.0 1.0,0.0 1.0 1.0,0.0 0.0 1.0)),\
             ((0.0 1.0 0.0,0.0 1.0 1.0,1.0 1.0 1.0,1.0 1.0 0.0,0.0 1.0 0.0))\
             )),\
             ((\
             ((0.0 0.0 0.0,0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 0.0,0.0 0.0 0.0)),\
             ((0.0 0.0 0.0,0.0 1.0 0.0,1.0 1.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0)),\
             ((0.0 0.0 0.0,1.0 0.0 0.0,1.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0)),\
             ((1.0 0.0 0.0,1.0 1.0 0.0,1.0 1.0 1.0,1.0 0.0 1.0,1.0 0.0 0.0)),\
             ((0.0 0.0 1.0,1.0 0.0 1.0,1.0 1.0 1.0,0.0 1.0 1.0,0.0 0.0 1.0)),\
             ((0.0 1.0 0.0,0.0 1.0 1.0,1.0 1.0 1.0,1.0 1.0 0.0,0.0 1.0 0.0))\
             ))\
             )",
        );
    }

    fn points_are_close2d(p1: (f64, f64), p2: (f64, f64)) -> bool {
        let tolerance = 1e-10;

        abs(p1.0 - p2.0) < tolerance && abs(p1.1 - p2.1) < tolerance
    }

    fn points_are_close3d(p1: (f64, f64, f64), p2: (f64, f64, f64)) -> bool {
        let tolerance = 1e-10;

        abs(p1.0 - p2.0) < tolerance && abs(p1.1 - p2.1) < tolerance && abs(p1.2 - p2.2) < tolerance
    }

    #[test]

    fn test_polygon_translation_2d() {
        let line = CoordSeq::<Point2d>::Linestring(vec![(0., 0.), (1.0, 1.0)])
            .to_sfcgal()
            .unwrap();

        let line = line.translate_2d(2., 3.).unwrap();

        let coords: CoordSeq<Point2d> = line.to_coordinates().unwrap();

        match coords {
            CoordSeq::Linestring(vec) => {
                assert!(points_are_close2d(vec[0], (2.0, 3.0)));

                assert!(points_are_close2d(vec[1], (3.0, 4.0)));
            }
            _ => panic!("Bad coordinates variant"),
        }
    }

    #[test]

    fn test_polygon_translation_3d() {
        let poly = CoordSeq::<Point3d>::Polygon(vec![vec![
            (0., 0., 0.),
            (1., 0., 0.),
            (1., 1., 0.),
            (0., 1., 0.),
            (0., 0., 0.),
        ]])
        .to_sfcgal()
        .unwrap();

        let poly = poly.translate_3d(1.0, 2.0, 3.0).unwrap();

        let coords: CoordSeq<Point3d> = poly.to_coordinates().unwrap();

        match coords {
            CoordSeq::Polygon(vec) => {
                let poly_test = &vec[0];

                assert!(points_are_close3d(poly_test[0], (1., 2., 3.)));

                assert!(points_are_close3d(poly_test[2], (2., 3., 3.)));
            }
            _ => panic!("Bad coordinate variant"),
        }
    }

    #[test]

    fn test_rotate_point() {
        let point = CoordSeq::Point((1., 0., 0.)).to_sfcgal().unwrap();

        let point = point.rotate_z(5. * PI / 2.).unwrap();

        let coords: CoordSeq<Point3d> = point.to_coordinates().unwrap();

        match coords {
            CoordSeq::Point(pt) => {
                assert!(points_are_close3d(pt, (0., 1., 0.)))
            }
            _ => panic!("Bad coordinate variant"),
        }
    }

    #[test]

    fn test_obj_export() {
        // ----------------------------------------
        // Open with Blender or any other 3d soft
        // to visually check
        // ----------------------------------------
        let temp_dir = env::temp_dir();

        let final_path = format!("{}/sfcgal_test.obj", temp_dir.to_str().unwrap());

        println!("Writing to {:?}", temp_dir);

        let input = "MULTISOLID Z (((((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),((1 1 0,0 1 0,0 1 1,1 1 1,1 1 0)),((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)))),((((2 4 6,2 5 6,3 5 6,3 4 6,2 4 6)),((2 4 7,3 4 7,3 5 7,2 5 7,2 4 7)),((2 4 6,3 4 6,3 4 7,2 4 7,2 4 6)),((3 5 6,2 5 6,2 5 7,3 5 7,3 5 6)),((3 4 6,3 5 6,3 5 7,3 4 7,3 4 6)),((2 4 6,2 4 7,2 5 7,2 5 6,2 4 6)))))";

        let shape = SFCGeometry::new(input).unwrap();

        assert!(shape.to_obj_file(final_path.as_str()).is_ok());
    }

    #[test]

    fn show_full_version() {
        println!("{:?}", _sfcgal_get_full_version());
    }
}
