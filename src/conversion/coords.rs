use crate::{
    utils::check_null_geom, GeomType, Point2d, Point3d, Result, SFCGeometry, ToCoordinates,
    ToSFCGAL,
};
use sfcgal_sys::{
    sfcgal_geometry_collection_add_geometry, sfcgal_geometry_collection_create,
    sfcgal_geometry_collection_geometry_n, sfcgal_geometry_collection_num_geometries,
    sfcgal_geometry_delete, sfcgal_geometry_t, sfcgal_linestring_add_point,
    sfcgal_linestring_create, sfcgal_linestring_num_points, sfcgal_linestring_point_n,
    sfcgal_multi_linestring_create, sfcgal_multi_point_create, sfcgal_multi_polygon_create,
    sfcgal_point_create_from_xy, sfcgal_point_create_from_xyz, sfcgal_point_x, sfcgal_point_y,
    sfcgal_point_z, sfcgal_polygon_add_interior_ring, sfcgal_polygon_create_from_exterior_ring,
    sfcgal_polygon_exterior_ring, sfcgal_polygon_interior_ring_n,
    sfcgal_polygon_num_interior_rings, sfcgal_polyhedral_surface_add_polygon,
    sfcgal_polyhedral_surface_create, sfcgal_polyhedral_surface_num_polygons,
    sfcgal_polyhedral_surface_polygon_n, sfcgal_solid_add_interior_shell, sfcgal_solid_create,
    sfcgal_solid_create_from_exterior_shell, sfcgal_solid_num_shells, sfcgal_solid_shell_n,
    sfcgal_triangle_create, sfcgal_triangle_set_vertex, sfcgal_triangle_vertex,
    sfcgal_triangulated_surface_add_triangle, sfcgal_triangulated_surface_create,
    sfcgal_triangulated_surface_num_triangles, sfcgal_triangulated_surface_triangle_n, size_t,
};

pub trait CoordType {}
impl CoordType for Point2d {}
impl CoordType for Point3d {}

/// Convert a point passed as a raw `sfcgal_geometry_t` to a tuple of coordinates.
pub trait FromSFCGALGeom {
    fn from_sfcgeometry(geom: *const sfcgal_geometry_t) -> Self;
}

/// Convert coordinates to `sfcgal_geometry_t`.
pub trait ToSFCGALGeom {
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t>;
}

impl FromSFCGALGeom for Point2d {
    fn from_sfcgeometry(c_geom: *const sfcgal_geometry_t) -> Point2d {
        let x = unsafe { sfcgal_point_x(c_geom) };
        let y = unsafe { sfcgal_point_y(c_geom) };
        (x, y)
    }
}

impl FromSFCGALGeom for Point3d {
    fn from_sfcgeometry(c_geom: *const sfcgal_geometry_t) -> Point3d {
        let x = unsafe { sfcgal_point_x(c_geom) };
        let y = unsafe { sfcgal_point_y(c_geom) };
        let z = unsafe { sfcgal_point_z(c_geom) };
        (x, y, z)
    }
}

impl ToSFCGALGeom for Point2d {
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t> {
        let g = unsafe { sfcgal_point_create_from_xy(self.0, self.1) };
        check_null_geom(g)?;
        Ok(g)
    }
}

impl ToSFCGALGeom for Point3d {
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t> {
        let g = unsafe { sfcgal_point_create_from_xyz(self.0, self.1, self.2) };
        check_null_geom(g)?;
        Ok(g)
    }
}

impl<T> ToSFCGALGeom for Vec<T>
where
    T: ToSFCGALGeom + CoordType,
{
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t> {
        let out_linestring = unsafe { sfcgal_linestring_create() };
        check_null_geom(out_linestring)?;
        for point in self {
            let sf_pt_geom: *mut sfcgal_geometry_t = point.to_sfcgeometry()?;
            unsafe { sfcgal_linestring_add_point(out_linestring, sf_pt_geom) };
        }
        Ok(out_linestring)
    }
}

/// Coordinates corresponding to the shapes described by SFCGAL Geometry types.
///
/// Used to construct SFCGeometry from coordinates or to retrieve coordinates from SFCGeometry.
/// ``` rust
/// use sfcgal::{CoordSeq, ToSFCGAL};
///
/// let coordinates = CoordSeq::Linestring(vec![(-0.5, -0.5, 2.5), (0., 0., 4.0)]);
/// let line_3d = coordinates.to_sfcgal().unwrap();
/// assert!(line_3d.is_valid().unwrap());
/// ```
/// ``` rust
/// use sfcgal::{SFCGeometry, CoordSeq, Point3d, ToCoordinates};
///
/// let line_3d = SFCGeometry::new("LINESTRING(3.5 5.6 1.0,4.8 10.5 1.0)").unwrap();
/// let coordinates = line_3d.to_coordinates::<Point3d>().unwrap();
/// // assert_eq!(coordinates, CoordSeq::Linestring(vec![(3.5, 5.6, 1.0), (4.8, 10.5, 1.0)]));
/// ```
#[derive(Debug, Clone, Hash, PartialEq)]
pub enum CoordSeq<T> {
    /// A Point is described by a tuple of 2 or 3 coordinates.
    Point(T),
    /// A Linestring is described by Vec of Point.
    Linestring(Vec<T>),
    /// A Polygon is described by a Vec of closed Linestring : the first one
    /// is the exterior ring and the others are interior rings.
    Polygon(Vec<Vec<T>>),
    /// A Multipoint is described by a Vec of Points.
    Multipoint(Vec<T>),
    /// A Multilinestring is described by a Vec of Linestrings.
    Multilinestring(Vec<Vec<T>>),
    /// A Multipolygon is described by a Vec of Polygons.
    Multipolygon(Vec<Vec<Vec<T>>>),
    /// A Geometrycollection is described by a Vec of any of these CoordSeq geometries.
    Geometrycollection(Vec<CoordSeq<T>>),
    /// A Polyhedralsurface is described by a Vec of Polygons.
    Polyhedralsurface(Vec<Vec<Vec<T>>>),
    /// A Triangulatedsurface is described by a Vec of Triangles.
    Triangulatedsurface(Vec<Vec<T>>),
    /// A Triangle is described by a Vec of 3 Points.
    Triangle(Vec<T>),
    /// A Solid is described by a Vec of Polyhedralsurfaces : the first one
    /// is the exterior shell and the others are interior shells.
    Solid(Vec<Vec<Vec<Vec<T>>>>),
    /// A Multisolid is described by a Vec of Solids.
    Multisolid(Vec<Vec<Vec<Vec<Vec<T>>>>>),
}

fn coords_polygon_to_sfcgal<T>(rings: &[Vec<T>]) -> Result<*mut sfcgal_geometry_t>
where
    T: ToSFCGALGeom + CoordType,
{
    let out_polygon =
        unsafe { sfcgal_polygon_create_from_exterior_ring(rings[0].to_sfcgeometry()?) };
    check_null_geom(out_polygon)?;
    for ring in rings.iter().skip(1) {
        unsafe { sfcgal_polygon_add_interior_ring(out_polygon, ring.to_sfcgeometry()?) };
    }
    Ok(out_polygon)
}

fn coords_polyhedralsurface_to_sfcgal<T>(
    polyhedres: &[Vec<Vec<T>>],
) -> Result<*mut sfcgal_geometry_t>
where
    T: ToSFCGALGeom + CoordType,
{
    let out_surf = unsafe { sfcgal_polyhedral_surface_create() };
    check_null_geom(out_surf)?;
    for poly in polyhedres {
        unsafe { sfcgal_polyhedral_surface_add_polygon(out_surf, coords_polygon_to_sfcgal(poly)?) };
    }
    Ok(out_surf)
}

unsafe fn pts_to_triangle_sfcgal<T>(pts: &[T]) -> Result<*mut sfcgal_geometry_t>
where
    T: ToSFCGALGeom,
{
    let out_triangle = sfcgal_triangle_create();
    check_null_geom(out_triangle)?;
    let p0: *mut sfcgal_geometry_t = pts[0].to_sfcgeometry()?;
    let p1: *mut sfcgal_geometry_t = pts[1].to_sfcgeometry()?;
    let p2: *mut sfcgal_geometry_t = pts[2].to_sfcgeometry()?;
    sfcgal_triangle_set_vertex(out_triangle, 0, p0);
    sfcgal_triangle_set_vertex(out_triangle, 1, p1);
    sfcgal_triangle_set_vertex(out_triangle, 2, p2);
    sfcgal_geometry_delete(p0);
    sfcgal_geometry_delete(p1);
    sfcgal_geometry_delete(p2);
    Ok(out_triangle)
}

/// Convert coordinates (tuple of 2 or 3 members) to [`SFCGeometry`] using
/// the corresponding [`CoordSeq`] variant of the wanted geometry.
///
/// [`SFCGeometry`]: struct.SFCGeometry.html
/// [`CoordSeq`]: enum.CoordSeq.html
impl<T: ToSFCGALGeom + CoordType> ToSFCGAL for CoordSeq<T> {
    /// Convert the coordinates of this [`CoordSeq`] to [`SFCGeometry`].
    ///
    /// [`SFCGeometry`]: struct.SFCGeometry.html
    /// [`CoordSeq`]: enum.CoordSeq.html
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        match self {
            CoordSeq::Point(pt) => unsafe { SFCGeometry::new_from_raw(pt.to_sfcgeometry()?, true) },
            CoordSeq::Multipoint(pts) => make_sfcgal_multi_geom!(
                sfcgal_multi_point_create(),
                pts.iter()
                    .map(|p| p.to_sfcgeometry())
                    .collect::<Result<Vec<_>>>()?
            ),
            CoordSeq::Linestring(pts) => {
                let out_linestring = pts.to_sfcgeometry()?;
                unsafe { SFCGeometry::new_from_raw(out_linestring, true) }
            }
            CoordSeq::Multilinestring(ref linestring_list) => make_sfcgal_multi_geom!(
                sfcgal_multi_linestring_create(),
                linestring_list
                    .iter()
                    .map(|l| l.to_sfcgeometry())
                    .collect::<Result<Vec<_>>>()?
            ),
            CoordSeq::Polygon(ref rings) => {
                let out_polygon = coords_polygon_to_sfcgal(rings)?;
                unsafe { SFCGeometry::new_from_raw(out_polygon, true) }
            }
            CoordSeq::Multipolygon(ref polygons) => make_sfcgal_multi_geom!(
                sfcgal_multi_polygon_create(),
                polygons
                    .iter()
                    .map(|p| coords_polygon_to_sfcgal(p))
                    .collect::<Result<Vec<_>>>()?
            ),
            CoordSeq::Geometrycollection(ref geoms) => {
                let out_geom_collection = unsafe { sfcgal_geometry_collection_create() };
                check_null_geom(out_geom_collection)?;
                for g_geom in geoms {
                    let mut sfcgeom = g_geom.to_sfcgal()?;
                    unsafe {
                        // ^^ sfcgeom is a SFCGeometry struct on which we own the inner C geometry,
                        // as we are passing it to `sfcgal_geometry_collection_add_geometry`
                        // we are not responsible for its deallocation anymore :
                        sfcgeom.owned = false;
                        sfcgal_geometry_collection_add_geometry(
                            out_geom_collection,
                            sfcgeom.c_geom.as_ptr(),
                        );
                    };
                }
                unsafe { SFCGeometry::new_from_raw(out_geom_collection, true) }
            }
            CoordSeq::Triangle(ref pts) => unsafe {
                SFCGeometry::new_from_raw(pts_to_triangle_sfcgal(pts)?, true)
            },
            CoordSeq::Triangulatedsurface(ref triangles) => unsafe {
                let out_surface = sfcgal_triangulated_surface_create();
                check_null_geom(out_surface)?;
                triangles
                    .iter()
                    .map(|pts| {
                        sfcgal_triangulated_surface_add_triangle(
                            out_surface,
                            pts_to_triangle_sfcgal(pts)?,
                        );
                        Ok(())
                    })
                    .collect::<Result<Vec<_>>>()?;
                SFCGeometry::new_from_raw(out_surface, true)
            },
            CoordSeq::Polyhedralsurface(ref polygons) => {
                let out_surface = coords_polyhedralsurface_to_sfcgal(polygons)?;
                unsafe { SFCGeometry::new_from_raw(out_surface, true) }
            }
            CoordSeq::Solid(ref polyhedres) => {
                let n_part = polyhedres.len();
                let out_solid = if n_part == 0 {
                    unsafe { sfcgal_solid_create() }
                } else if n_part == 1 {
                    let exterior = coords_polyhedralsurface_to_sfcgal(&polyhedres[0])?;
                    unsafe { sfcgal_solid_create_from_exterior_shell(exterior) }
                } else {
                    let exterior = coords_polyhedralsurface_to_sfcgal(&polyhedres[0])?;
                    let r_solid = unsafe { sfcgal_solid_create_from_exterior_shell(exterior) };
                    polyhedres
                        .into_iter()
                        .skip(1)
                        .map(|poly| {
                            unsafe {
                                sfcgal_solid_add_interior_shell(
                                    r_solid,
                                    coords_polyhedralsurface_to_sfcgal(poly)?,
                                )
                            };
                            Ok(())
                        })
                        .collect::<Result<Vec<_>>>()?;
                    r_solid
                };
                unsafe { SFCGeometry::new_from_raw(out_solid, true) }
            }
            _ => unimplemented!(),
        }
    }
}

macro_rules! sfcgal_pt_to_coords {
    ($geom: expr, T) => {{
        T::from_sfcgeometry($geom)
    }};
}

macro_rules! sfcgal_line_to_coords {
    ($geom: expr, $type_pt: ident) => {{
        let g = $geom;
        let n_points = unsafe { sfcgal_linestring_num_points(g) };
        let mut v_points = Vec::with_capacity(n_points as usize);
        for i in 0..n_points {
            let pt_sfcgal = unsafe { sfcgal_linestring_point_n(g, i) };
            // check_null_geom(g)?;
            v_points.push(sfcgal_pt_to_coords!(pt_sfcgal, $type_pt));
        }
        v_points
    }};
}

macro_rules! sfcgal_polygon_to_coords {
    ($geom: expr, $type_pt: ident) => {{
        let g = $geom;
        let nrings = unsafe { sfcgal_polygon_num_interior_rings(g) };
        let exterior_sfcgal = unsafe { sfcgal_polygon_exterior_ring(g) };
        let mut rings = Vec::with_capacity(nrings as usize + 1);
        rings.push(sfcgal_line_to_coords!(exterior_sfcgal, $type_pt));
        for ix in 0..nrings {
            let line_sfcgal = unsafe { sfcgal_polygon_interior_ring_n(g, ix) };
            rings.push(sfcgal_line_to_coords!(line_sfcgal, $type_pt));
        }
        rings
    }};
}

macro_rules! sfcgal_triangle_to_coords {
    ($geom: expr, $type_pt: ident) => {{
        let g = $geom;
        let (p0, p1, p2) = unsafe {
            (
                sfcgal_pt_to_coords!(sfcgal_triangle_vertex(g, 0), $type_pt),
                sfcgal_pt_to_coords!(sfcgal_triangle_vertex(g, 1), $type_pt),
                sfcgal_pt_to_coords!(sfcgal_triangle_vertex(g, 2), $type_pt),
            )
        };
        vec![p0, p1, p2]
    }};
}

macro_rules! sfcgal_polyhedral_surface_to_coords {
    ($geom: expr, $type_pt: ident) => {{
        let g = $geom;
        let ngeoms = unsafe { sfcgal_polyhedral_surface_num_polygons(g) };
        let mut polygons = Vec::with_capacity(ngeoms as usize);
        for ix in 0..ngeoms {
            let poly = unsafe { sfcgal_polyhedral_surface_polygon_n(g, ix) };
            polygons.push(sfcgal_polygon_to_coords!(poly, $type_pt));
        }
        polygons
    }};
}

fn get_nb_geometry(geom: *const sfcgal_geometry_t) -> usize {
    unsafe { sfcgal_geometry_collection_num_geometries(geom) as usize }
}

fn get_geom_at_index(geom: *const sfcgal_geometry_t, ix: usize) -> *const sfcgal_geometry_t {
    unsafe { sfcgal_geometry_collection_geometry_n(geom, ix as size_t) }
}

/// Convert a [`SFCGeometry`], given it's internal [`GeomType`], to the corresponding [`CoordSeq`]
/// holding its coordinates (as tuple of 2 or 3 members).
///
/// [`SFCGeometry`]: struct.SFCGeometry.html
/// [`GeomType`]: enum.GeomType.html
/// [`CoordSeq`]: enum.CoordSeq.html
impl ToCoordinates for SFCGeometry {
    fn to_coordinates<T>(&self) -> Result<CoordSeq<T>>
    where
        T: FromSFCGALGeom + CoordType,
    {
        match self._type()? {
            GeomType::Point => Ok(CoordSeq::Point(sfcgal_pt_to_coords!(
                unsafe { self.c_geom.as_ref() },
                T
            ))),
            GeomType::Multipoint => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = get_nb_geometry(geom);
                let mut pts = Vec::with_capacity(ngeoms);
                for ix in 0..ngeoms {
                    pts.push(sfcgal_pt_to_coords!(get_geom_at_index(geom, ix), T));
                }
                Ok(CoordSeq::Multipoint(pts))
            }
            GeomType::Linestring => Ok(CoordSeq::Linestring(sfcgal_line_to_coords!(
                unsafe { self.c_geom.as_ref() },
                T
            ))),
            GeomType::Multilinestring => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = get_nb_geometry(geom);
                let mut lines = Vec::with_capacity(ngeoms);
                for ix in 0..ngeoms {
                    lines.push(sfcgal_line_to_coords!(get_geom_at_index(geom, ix), T));
                }
                Ok(CoordSeq::Multilinestring(lines))
            }
            GeomType::Polygon => {
                let geom = unsafe { self.c_geom.as_ref() };
                Ok(CoordSeq::Polygon(sfcgal_polygon_to_coords!(geom, T)))
            }
            GeomType::Multipolygon => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = get_nb_geometry(geom);
                let mut polygons = Vec::with_capacity(ngeoms);
                for ix in 0..ngeoms {
                    polygons.push(sfcgal_polygon_to_coords!(get_geom_at_index(geom, ix), T));
                }
                Ok(CoordSeq::Multipolygon(polygons))
            }
            GeomType::Geometrycollection => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = get_nb_geometry(geom);
                let mut geoms = Vec::with_capacity(ngeoms);
                for ix in 0..ngeoms {
                    let inner_geom = unsafe {
                        // Todo : document what we are doing that
                        SFCGeometry::new_from_raw(
                            sfcgal_geometry_collection_geometry_n(geom, ix as size_t)
                                as *mut sfcgal_geometry_t,
                            false,
                        )?
                    };
                    let coords: CoordSeq<T> = inner_geom.to_coordinates()?;
                    geoms.push(coords);
                }
                Ok(CoordSeq::Geometrycollection(geoms))
            }
            GeomType::Triangle => Ok(CoordSeq::Triangle(sfcgal_triangle_to_coords!(
                unsafe { self.c_geom.as_ref() },
                T
            ))),
            GeomType::Triangulatedsurface => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = unsafe { sfcgal_triangulated_surface_num_triangles(geom) };
                let mut triangles = Vec::with_capacity(ngeoms as usize);
                for ix in 0..ngeoms {
                    let triangle = unsafe { sfcgal_triangulated_surface_triangle_n(geom, ix) };
                    triangles.push(sfcgal_triangle_to_coords!(triangle, T));
                }
                Ok(CoordSeq::Triangulatedsurface(triangles))
            }
            GeomType::Polyhedralsurface => Ok(CoordSeq::Polyhedralsurface(
                sfcgal_polyhedral_surface_to_coords!(unsafe { self.c_geom.as_ref() }, T),
            )),
            GeomType::Solid => {
                let geom = unsafe { self.c_geom.as_ref() };
                let ngeoms = unsafe { sfcgal_solid_num_shells(geom) };
                let mut polyhedres = Vec::with_capacity(ngeoms as usize);
                for ix in 0..ngeoms {
                    let poly = unsafe { sfcgal_solid_shell_n(geom, ix) };
                    polyhedres.push(sfcgal_polyhedral_surface_to_coords!(poly, T));
                }
                Ok(CoordSeq::Solid(polyhedres))
            }
            GeomType::Multisolid => {
                let geom_ms = unsafe { self.c_geom.as_ref() };
                let n_solides = get_nb_geometry(geom_ms);
                let mut solides = Vec::with_capacity(n_solides);
                for ix_geom in 0..n_solides {
                    let solid = get_geom_at_index(geom_ms, ix_geom);
                    let n_shell = unsafe { sfcgal_solid_num_shells(solid) };
                    let mut polyhedres = Vec::with_capacity(n_shell as usize);
                    for ix in 0..n_shell {
                        let poly = unsafe { sfcgal_solid_shell_n(solid, ix) };
                        polyhedres.push(sfcgal_polyhedral_surface_to_coords!(poly, T));
                    }
                    solides.push(polyhedres);
                }
                Ok(CoordSeq::Multisolid(solides))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Point2d, Point3d};
    use crate::*;

    #[test]
    fn point_2d_sfcgal_to_coordinates_to_sfcgal() {
        let input_wkt = "POINT(0.1 0.9)";
        let pt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point2d> = pt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Point(ref coords) = coords {
            assert_ulps_eq!(coords.0, 0.1);
            assert_ulps_eq!(coords.1, 0.9);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let pt_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, pt_sfcgal.to_wkt_decim(1).unwrap())
    }

    #[test]
    fn point_3d_sfcgal_to_coordinates_to_sfcgal() {
        let input_wkt = "POINT(0.1 0.9 1.0)";
        let pt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = pt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Point(ref coords) = coords {
            assert_ulps_eq!(coords.0, 0.1);
            assert_ulps_eq!(coords.1, 0.9);
            assert_ulps_eq!(coords.2, 1.0);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let pt_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, pt_sfcgal.to_wkt_decim(1).unwrap())
    }

    #[test]
    fn multipoint_2d_sfcgal_to_coordinates_to_sfcgal() {
        let input_wkt = "MULTIPOINT((3.5 5.6),(4.8 10.5))";
        let multipt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point2d> = multipt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multipoint(ref coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let mpt_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, mpt_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn multipoint_3d_sfcgal_to_coordinates_to_sfcgal() {
        let input_wkt = "MULTIPOINT((3.5 5.6 1.0),(4.8 10.5 1.0))";
        let multipt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = multipt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multipoint(ref coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[0].2, 1.0);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
            assert_ulps_eq!(coords[1].2, 1.0);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let mpt_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, mpt_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn linestring_2d_sfcgal_to_coordinates() {
        let input_wkt = "LINESTRING(3.5 5.6,4.8 10.5)";
        let linestring_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point2d> = linestring_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Linestring(ref coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let line_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, line_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn linestring_3d_sfcgal_to_coordinates() {
        let input_wkt = "LINESTRING(3.5 5.6 1.0,4.8 10.5 1.0)";
        let linestring_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = linestring_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Linestring(ref coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[0].2, 1.0);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
            assert_ulps_eq!(coords[1].2, 1.0);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let line_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, line_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn multilinestring_2d_sfcgal_to_coordinates() {
        let input_wkt = "MULTILINESTRING((3.5 5.6,4.8 10.5),(8.9 12.9,2.1 3.5),(1.1 4.8,6.2 7.5))";
        let mls_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point2d> = mls_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multilinestring(ref coords) = coords {
            assert_ulps_eq!(coords[0][0].0, 3.5);
            assert_ulps_eq!(coords[0][0].1, 5.6);
            assert_ulps_eq!(coords[1][1].0, 2.1);
            assert_ulps_eq!(coords[1][1].1, 3.5);
            assert_ulps_eq!(coords[2][0].0, 1.1);
            assert_ulps_eq!(coords[2][0].1, 4.8);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let mls_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, mls_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn multilinestring_3d_sfcgal_to_coordinates() {
        let input_wkt = "MULTILINESTRING((3.5 5.6 1.0,4.8 10.5 1.0),(8.9 12.9 4.0,2.1 3.5 4.0),(1.1 4.8 4.0,6.2 7.5 4.0))";
        let multilinestring_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = multilinestring_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multilinestring(ref coords) = coords {
            assert_ulps_eq!(coords[0][0].0, 3.5);
            assert_ulps_eq!(coords[0][0].1, 5.6);
            assert_ulps_eq!(coords[0][0].2, 1.0);
            assert_ulps_eq!(coords[1][1].0, 2.1);
            assert_ulps_eq!(coords[1][1].1, 3.5);
            assert_ulps_eq!(coords[1][1].2, 4.0);
            assert_ulps_eq!(coords[2][0].0, 1.1);
            assert_ulps_eq!(coords[2][0].1, 4.8);
            assert_ulps_eq!(coords[2][0].2, 4.0);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let mls_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, mls_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn triangle_2d_sfcgal_to_coordinates_to_sfcgal() {
        let input_wkt = "TRIANGLE((0.0 1.0,1.0 0.0,1.0 1.0,0.0 1.0))";
        let triangle_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point2d> = triangle_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Triangle(ref pts) = coords {
            assert_ulps_eq!(pts[0].0, 0.0);
            assert_ulps_eq!(pts[0].1, 1.0);
            assert_ulps_eq!(pts[1].0, 1.0);
            assert_ulps_eq!(pts[1].1, 0.0);
            assert_ulps_eq!(pts[2].0, 1.0);
            assert_ulps_eq!(pts[2].1, 1.0);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let triangle_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, triangle_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn solid_3d_sfcgal_to_coordinates() {
        let input_wkt = "SOLID((((0 0 0,0 0 1,0 1 1,0 1 0,0 0 0)),\
                         ((0 0 0,0 1 0,1 1 0,1 0 0,0 0 0)),\
                         ((0 0 0,1 0 0,1 0 1,0 0 1,0 0 0)),\
                         ((1 0 0,1 1 0,1 1 1,1 0 1,1 0 0)),\
                         ((0 0 1,1 0 1,1 1 1,0 1 1,0 0 1)),\
                         ((0 1 0,0 1 1,1 1 1,1 1 0,0 1 0))))";
        let cube = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = cube.to_coordinates().unwrap();
        if let CoordSeq::Solid(ref polys) = coords {
            let coords_1st_face_polygon = &polys[0][0];
            assert_eq!(
                format!("{:?}", coords_1st_face_polygon),
                "[[(0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)]]");
            assert_eq!(polys[0].len(), 6usize);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let cube_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, cube_sfcgal.to_wkt_decim(0).unwrap());
    }

    #[test]
    fn solid_with_interior_shell_3d_sfcgal_to_coordinates() {
        let input_wkt = "SOLID((\
                         ((0.0 0.0 0.0,0.0 1.0 0.0,1.0 1.0 0.0,1.0 0.0 0.0,0.0 0.0 0.0)),\
                         ((1.0 0.0 0.0,1.0 1.0 0.0,1.0 1.0 1.0,1.0 0.0 1.0,1.0 0.0 0.0)),\
                         ((0.0 1.0 0.0,0.0 1.0 1.0,1.0 1.0 1.0,1.0 1.0 0.0,0.0 1.0 0.0)),\
                         ((0.0 0.0 1.0,0.0 1.0 1.0,0.0 1.0 0.0,0.0 0.0 0.0,0.0 0.0 1.0)),\
                         ((1.0 0.0 1.0,1.0 1.0 1.0,0.0 1.0 1.0,0.0 0.0 1.0,1.0 0.0 1.0)),\
                         ((1.0 0.0 0.0,1.0 0.0 1.0,0.0 0.0 1.0,0.0 0.0 0.0,1.0 0.0 0.0))\
                         ),(\
                         ((0.0 0.0 0.0,0.0 0.5 0.0,0.5 0.5 0.0,0.5 0.0 0.0,0.0 0.0 0.0)),\
                         ((0.5 0.0 0.0,0.5 0.5 0.0,0.5 0.5 0.5,0.5 0.0 0.5,0.5 0.0 0.0)),\
                         ((0.0 0.5 0.0,0.0 0.5 0.5,0.5 0.5 0.5,0.5 0.5 0.0,0.0 0.5 0.0)),\
                         ((0.0 0.0 0.5,0.0 0.5 0.5,0.0 0.5 0.0,0.0 0.0 0.0,0.0 0.0 0.5)),\
                         ((0.5 0.0 0.5,0.5 0.5 0.5,0.0 0.5 0.5,0.0 0.0 0.5,0.5 0.0 0.5)),\
                         ((0.5 0.0 0.0,0.5 0.0 0.5,0.0 0.0 0.5,0.0 0.0 0.0,0.5 0.0 0.0))\
                         ))";
        let cube = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = cube.to_coordinates().unwrap();
        if let CoordSeq::Solid(ref polys) = coords {
            let num_shell = polys.len();
            assert_eq!(num_shell, 2usize);
            let coords_1st_face_polygon = &polys[0][0];
            assert_eq!(
                format!("{:?}", coords_1st_face_polygon),
                "[[(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 0.0)]]");
            assert_eq!(polys[0].len(), 6usize);
        } else {
            panic!("Unexpected coordinates when converting from SFCGAL to coordinates as tuples.")
        }
        let cube_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, cube_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn geometrycollection_3d_sfcgal_to_coordinates() {
        let input_wkt =
            "GEOMETRYCOLLECTION(POINT(4.0 6.0 1.0),LINESTRING(4.0 6.0 1.0,7.0 10.0 1.0))";
        let gc_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let coords: CoordSeq<Point3d> = gc_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Geometrycollection(ref seq_coords) = coords {
            if let CoordSeq::Point(ref coords) = &seq_coords[0] {
                assert_ulps_eq!(coords.0, 4.0);
                assert_ulps_eq!(coords.1, 6.0);
                assert_ulps_eq!(coords.2, 1.0);
            } else {
                panic!("Error with first member of geometry collection");
            }
            if let CoordSeq::Linestring(ref coords) = &seq_coords[1] {
                assert_ulps_eq!(coords[0].0, 4.0);
                assert_ulps_eq!(coords[0].1, 6.0);
                assert_ulps_eq!(coords[0].2, 1.0);
                assert_ulps_eq!(coords[1].0, 7.0);
                assert_ulps_eq!(coords[1].1, 10.0);
                assert_ulps_eq!(coords[1].2, 1.0);
            } else {
                panic!("Error with second member of geometry collection");
            }
        }
        let gc_sfcgal = coords.to_sfcgal().unwrap();
        assert_eq!(input_wkt, gc_sfcgal.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn sfcgal_geometry_construction_from_coordinates() {
        let coords_point = (1.1, 1.1);
        let coords_multipoint = vec![(1.1, 1.1), (2.2, 2.2)];
        let coords_linestring = vec![(1.1, 1.1), (2.2, 2.2)];
        let coords_multilinestring = vec![
            vec![(1.1, 1.1), (2.2, 2.2)],
            vec![(3.1, 3.1), (5.2, 5.2), (5.2, 5.2)],
            vec![(1.1, 1.1), (2.2, 2.2), (5.2, 5.2)],
        ];
        let coords_polygon = vec![
            vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.), (0., 0.)], // Exterior ring
            vec![(0.1, 0.1), (0.1, 0.9), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)], // 1 interior ring
        ];
        let coords_multipolygon = vec![
            vec![
                // First polygon
                vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.), (0., 0.)], // Exterior ring
                vec![(0.1, 0.1), (0.1, 0.9), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)], // 1 interior ring
            ],
            vec![
                // second polygon
                vec![(2., 2.), (6., 2.), (6., 6.), (2., 6.), (2., 2.)], // Exterior ring
            ],
        ];
        let pt_sfcgal = CoordSeq::Point(coords_point).to_sfcgal().unwrap();
        let multipt_sfcgal = CoordSeq::Multipoint(coords_multipoint).to_sfcgal().unwrap();
        let line_sfcgal = CoordSeq::Linestring(coords_linestring).to_sfcgal().unwrap();
        let multiline_sfcgal = CoordSeq::Multilinestring(coords_multilinestring)
            .to_sfcgal()
            .unwrap();
        let polygon_sfcgal = CoordSeq::Polygon(coords_polygon).to_sfcgal().unwrap();
        let multipolygon_sfcgal = CoordSeq::Multipolygon(coords_multipolygon)
            .to_sfcgal()
            .unwrap();
        assert_eq!("POINT(1.1 1.1)", pt_sfcgal.to_wkt_decim(1).unwrap());
        assert_eq!(
            "MULTIPOINT((1.1 1.1),(2.2 2.2))",
            multipt_sfcgal.to_wkt_decim(1).unwrap()
        );
        assert_eq!(
            "LINESTRING(1.1 1.1,2.2 2.2)",
            line_sfcgal.to_wkt_decim(1).unwrap()
        );
        assert_eq!(
            "MULTILINESTRING((1.1 1.1,2.2 2.2),\
             (3.1 3.1,5.2 5.2,5.2 5.2),\
             (1.1 1.1,2.2 2.2,5.2 5.2))",
            multiline_sfcgal.to_wkt_decim(1).unwrap(),
        );
        assert_eq!(
            "POLYGON((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0),\
             (0.1 0.1,0.1 0.9,0.9 0.9,0.9 0.1,0.1 0.1))",
            polygon_sfcgal.to_wkt_decim(1).unwrap(),
        );
        assert_eq!(
            "MULTIPOLYGON(((0.0 0.0,1.0 0.0,1.0 1.0,0.0 1.0,0.0 0.0),\
             (0.1 0.1,0.1 0.9,0.9 0.9,0.9 0.1,0.1 0.1)),\
             ((2.0 2.0,6.0 2.0,6.0 6.0,2.0 6.0,2.0 2.0)))",
            multipolygon_sfcgal.to_wkt_decim(1).unwrap(),
        );
    }

    #[test]
    fn sfcgal_geometry_construction_from_coordinates_3d() {
        let coords_point = (1.1, 1.1, 1.1);
        let coords_multipoint = vec![(1.1, 1.1, 5.0), (2.2, 2.2, 5.0)];
        let coords_linestring = vec![(1.1, 1.1, 5.0), (2.2, 2.2, 5.0)];
        let coords_multilinestring = vec![
            vec![(1.1, 1.1, 3.0), (2.2, 2.2, 3.0)],
            vec![(3.1, 3.1, 3.0), (5.2, 5.2, 3.0), (5.2, 5.2, 3.0)],
            vec![(1.1, 1.1, 3.0), (2.2, 2.2, 3.0), (5.2, 5.2, 3.0)],
        ];
        let coords_polygon = vec![
            vec![
                (0., 0., 3.0),
                (1., 0., 3.0),
                (1., 1., 3.0),
                (0., 1., 3.0),
                (0., 0., 3.0),
            ], // Exterior ring
            vec![
                (0.1, 0.1, 3.0),
                (0.1, 0.9, 3.0),
                (0.9, 0.9, 3.0),
                (0.9, 0.1, 3.0),
                (0.1, 0.1, 3.0),
            ], // 1 interior ring
        ];
        let coords_multipolygon = vec![
            vec![
                // First polygon
                vec![
                    (0., 0., 3.0),
                    (1., 0., 3.0),
                    (1., 1., 3.0),
                    (0., 1., 3.0),
                    (0., 0., 3.0),
                ], // Exterior ring
                vec![
                    (0.1, 0.1, 3.0),
                    (0.1, 0.9, 3.0),
                    (0.9, 0.9, 3.0),
                    (0.9, 0.1, 3.0),
                    (0.1, 0.1, 3.0),
                ], // 1 interior ring
            ],
            vec![
                // second polygon
                vec![
                    (2., 2., 3.),
                    (6., 2., 3.),
                    (6., 6., 3.),
                    (2., 6., 3.),
                    (2., 2., 3.),
                ], // Exterior ring
            ],
        ];
        let pt_sfcgal = CoordSeq::Point(coords_point).to_sfcgal().unwrap();
        let multipt_sfcgal = CoordSeq::Multipoint(coords_multipoint).to_sfcgal().unwrap();
        let line_sfcgal = CoordSeq::Linestring(coords_linestring).to_sfcgal().unwrap();
        let multiline_sfcgal = CoordSeq::Multilinestring(coords_multilinestring)
            .to_sfcgal()
            .unwrap();
        let polygon_sfcgal = CoordSeq::Polygon(coords_polygon).to_sfcgal().unwrap();
        let multipolygon_sfcgal = CoordSeq::Multipolygon(coords_multipolygon)
            .to_sfcgal()
            .unwrap();
        assert_eq!("POINT(1.1 1.1 1.1)", pt_sfcgal.to_wkt_decim(1).unwrap());
        assert_eq!(
            "MULTIPOINT((1.1 1.1 5.0),(2.2 2.2 5.0))",
            multipt_sfcgal.to_wkt_decim(1).unwrap()
        );
        assert_eq!(
            "LINESTRING(1.1 1.1 5.0,2.2 2.2 5.0)",
            line_sfcgal.to_wkt_decim(1).unwrap()
        );
        assert_eq!(
            "MULTILINESTRING((1.1 1.1 3.0,2.2 2.2 3.0),(3.1 3.1 3.0,5.2 5.2 3.0,5.2 5.2 3.0),\
             (1.1 1.1 3.0,2.2 2.2 3.0,5.2 5.2 3.0))",
            multiline_sfcgal.to_wkt_decim(1).unwrap(),
        );
        assert_eq!(
            "POLYGON((0.0 0.0 3.0,1.0 0.0 3.0,1.0 1.0 3.0,0.0 1.0 3.0,0.0 0.0 3.0),\
             (0.1 0.1 3.0,0.1 0.9 3.0,0.9 0.9 3.0,0.9 0.1 3.0,0.1 0.1 3.0))",
            polygon_sfcgal.to_wkt_decim(1).unwrap(),
        );
        assert_eq!(
            "MULTIPOLYGON(((0.0 0.0 3.0,1.0 0.0 3.0,1.0 1.0 3.0,0.0 1.0 3.0,0.0 0.0 3.0),\
             (0.1 0.1 3.0,0.1 0.9 3.0,0.9 0.9 3.0,0.9 0.1 3.0,0.1 0.1 3.0)),\
             ((2.0 2.0 3.0,6.0 2.0 3.0,6.0 6.0 3.0,2.0 6.0 3.0,2.0 2.0 3.0)))",
            multipolygon_sfcgal.to_wkt_decim(1).unwrap(),
        );
    }

}
