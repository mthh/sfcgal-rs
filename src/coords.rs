use failure::Error;
#[allow(unused_imports)]
use sfcgal_sys::{
    sfcgal_geometry_t,
    sfcgal_point_create_from_xy, sfcgal_point_x, sfcgal_point_y, sfcgal_point_z,
    sfcgal_linestring_create, sfcgal_linestring_add_point, sfcgal_linestring_point_n, sfcgal_linestring_num_points,
    sfcgal_triangle_create, sfcgal_triangle_set_vertex_from_xy, sfcgal_triangle_vertex,
    sfcgal_polygon_create, sfcgal_polygon_create_from_exterior_ring, sfcgal_polygon_add_interior_ring,
    sfcgal_polygon_num_interior_rings, sfcgal_polygon_exterior_ring, sfcgal_polygon_interior_ring_n,
    sfcgal_multi_point_create, sfcgal_multi_linestring_create, sfcgal_multi_polygon_create,
    sfcgal_geometry_collection_add_geometry, sfcgal_geometry_collection_num_geometries,
    sfcgal_geometry_collection_geometry_n, sfcgal_geometry_collection_create,
};
use crate::{SFCGeometry, GeomType, Result, ToCoordinates};

/// Type for manipulating 2d coordinates, reprensented as (x, y).
pub type point_2d = (f64, f64);

/// Type for manipulating 3d coordinates, reprensented as (x, y, z).
pub type point_3d = (f64, f64, f64);

/// Convert a point passed as a raw `sfcgal_geometry_t` to a tuple of coordinates.
pub trait FromSfcgalGeom {
    fn from_sfcgeometry(geom: *const sfcgal_geometry_t) -> Self;
}

impl FromSfcgalGeom for point_2d {
    fn from_sfcgeometry(c_geom: *const sfcgal_geometry_t) -> point_2d {
        let x = unsafe { sfcgal_point_x(c_geom) };
        let y = unsafe { sfcgal_point_y(c_geom) };
        (x, y)
    }
}

impl FromSfcgalGeom for point_3d {
    fn from_sfcgeometry(c_geom: *const sfcgal_geometry_t) -> point_3d {
        let x = unsafe { sfcgal_point_x(c_geom) };
        let y = unsafe { sfcgal_point_y(c_geom) };
        let z = unsafe { sfcgal_point_z(c_geom) };
        (x, y, z)
    }
}
/// List of coordinates corresponding to the shapes described by SFCGAL Geometry types.
#[derive(Debug, Clone, Hash)]
pub enum CoordSeq<T> {
    Point(T),
    Linestring(Vec<T>),
    Polygon(Vec<Vec<T>>),
    Multipoint(Vec<T>),
    Multilinestring(Vec<Vec<T>>),
    Multipolygon(Vec<Vec<Vec<T>>>),
    Geometrycollection(Vec<CoordSeq<T>>),
    Polyhedralsurface(Vec<Vec<Vec<T>>>),
    Triangulatedsurface(Vec<Vec<T>>),
    Triangle(Vec<T>),
    Solid(Vec<Vec<Vec<T>>>),
    Multisolid(Vec<Vec<Vec<Vec<T>>>>),
}

macro_rules! sfcgal_pt_to_coord {
    ($geom: expr, point_2d) => ({
        point_2d::from_sfcgeometry($geom)
    });
    ($geom: expr, point_3d) => ({
        point_3d::from_sfcgeometry($geom)
    });
    ($geom: expr, T) => ({
        T::from_sfcgeometry($geom)
    });
}

macro_rules! sfcgal_line_to_coord {
    ($geom: expr, $type_pt: ident) => ({
        let g = $geom;
        let n_points = unsafe { sfcgal_linestring_num_points(g) };
        let mut v_points = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let pt_sfcgal = unsafe { sfcgal_linestring_point_n(g, i) };
            // check_null_geom(g)?;
            v_points.push(sfcgal_pt_to_coord!(pt_sfcgal, $type_pt));
        }
        v_points
    });
}

macro_rules! get_nb_geometry {
    ($geom: expr) => ({
        unsafe {
            sfcgal_geometry_collection_num_geometries($geom)
        }
    });
}

macro_rules! get_geom_at_index {
    ($geom: expr, $ix: expr) => ({
        unsafe {
            sfcgal_geometry_collection_geometry_n($geom, $ix)
        }
    });
}


impl<T: FromSfcgalGeom> ToCoordinates<T> for SFCGeometry {
    fn to_coordinates(&self) ->  Result<CoordSeq<T>> {
        match self._type()? {
            GeomType::Point => {
                Ok(
                    CoordSeq::Point(
                        sfcgal_pt_to_coord!(unsafe { self.c_geom.as_ref() }, T)
                    )
                )
            },
            GeomType::Multipoint => {
                let ngeoms = get_nb_geometry!(self.c_geom.as_ref());
                let mut pts = Vec::with_capacity(ngeoms);
                for ix in 0..ngeoms {
                    pts.push(
                        sfcgal_pt_to_coord!(get_geom_at_index!(self.c_geom.as_ref(), ix), T));
                }
                Ok(CoordSeq::Multipoint(pts))
            },
            GeomType::Linestring => {
                Ok(
                    CoordSeq::Linestring(
                        sfcgal_line_to_coord!(unsafe { self.c_geom.as_ref() }, T)
                    )
                )
            },
            // GeomType::Multilinestring => { },
            // GeomType::Polygon => { },
            // GeomType::Multipolygon => { },
            // GeomType::Geometrycollection => { },
            // GeomType::Polyhedralsurface => { },
            // GeomType::Triangulatedsurface => { },
            // GeomType::Triangle => { },
            // GeomType::Solid => { },
            // GeomType::Multisolid => { },
            _ => unimplemented!()
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use super::{point_2d, point_3d};

    #[test]
    fn point_2d_sfcgal_to_coordinates() {
        let pt_sfcgal = SFCGeometry::new("POINT(0.1 0.9)").unwrap();
        let coords: CoordSeq<point_2d> = pt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Point(coords) = coords {
            assert_ulps_eq!(coords.0, 0.1);
            assert_ulps_eq!(coords.1, 0.9);
        }
    }

    #[test]
    fn point_3d_sfcgal_to_coordinates() {
        let pt_sfcgal = SFCGeometry::new("POINT(0.1 0.9 1.0)").unwrap();
        let coords: CoordSeq<point_3d> = pt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Point(coords) = coords {
            assert_ulps_eq!(coords.0, 0.1);
            assert_ulps_eq!(coords.1, 0.9);
            assert_ulps_eq!(coords.2, 1.0);
        }
    }

    #[test]
    fn multipoint_2d_sfcgal_to_coordinates() {
        let multipt_sfcgal = SFCGeometry::new("MULTIPOINT((3.5 5.6), (4.8 10.5))").unwrap();
        let coords: CoordSeq<point_2d> = multipt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multipoint(coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
        }
    }

    #[test]
    fn multipoint_3d_sfcgal_to_coordinates() {
        let multipt_sfcgal = SFCGeometry::new("MULTIPOINT((3.5 5.6 1.0), (4.8 10.5 1.0))").unwrap();
        let coords: CoordSeq<point_3d> = multipt_sfcgal.to_coordinates().unwrap();
        if let CoordSeq::Multipoint(coords) = coords {
            assert_ulps_eq!(coords[0].0, 3.5);
            assert_ulps_eq!(coords[0].1, 5.6);
            assert_ulps_eq!(coords[0].2, 1.0);
            assert_ulps_eq!(coords[1].0, 4.8);
            assert_ulps_eq!(coords[1].1, 10.5);
            assert_ulps_eq!(coords[1].2, 1.0);
        }
    }
}
