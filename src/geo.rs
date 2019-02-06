use failure::Error;
use sfcgal_sys::{
    sfcgal_geometry_t,
    sfcgal_point_create_from_xy, sfcgal_point_x, sfcgal_point_y, sfcgal_point_z,
    sfcgal_linestring_create, sfcgal_linestring_add_point, sfcgal_linestring_point_n, sfcgal_linestring_num_points,
    sfcgal_polygon_create, sfcgal_polygon_create_from_exterior_ring, sfcgal_polygon_add_interior_ring,
    sfcgal_polygon_num_interior_rings, sfcgal_polygon_exterior_ring, sfcgal_polygon_interior_ring_n,
    sfcgal_multi_point_create, sfcgal_multi_linestring_create, sfcgal_multi_polygon_create,
    sfcgal_geometry_collection_add_geometry, sfcgal_geometry_collection_num_geometries,
    sfcgal_geometry_collection_geometry_n, sfcgal_geometry_collection_create,
};
use crate::{Result, SFCGeometry, GeomType, ToSfcgal, utils::check_null_geom};

// define our own TryInto while the std trait is not stable
pub trait TryInto<T> {
    type Err;
    fn try_into(self) -> Result<T>;
}

impl<'a> TryInto<geo_types::Geometry<f64>> for SFCGeometry {
    type Err = Error;

    fn try_into(self) -> Result<geo_types::Geometry<f64>> {
        match self._type()? {
            GeomType::Point => {
                Ok(
                    geo_types::Geometry::Point(
                        unsafe { geo_point_from_sfcgal(self.0.as_ref()) }
                    )
                )
            },
            GeomType::Multipoint => {
                let ngeoms = unsafe {
                    sfcgal_geometry_collection_num_geometries(self.0.as_ref())
                };
                let mut pts = Vec::with_capacity(ngeoms);
                for i in 0..ngeoms {
                    let geom = unsafe { sfcgal_geometry_collection_geometry_n(self.0.as_ref(), i) };
                    pts.push(geo_point_from_sfcgal(geom));
                }
                Ok(
                    geo_types::Geometry::MultiPoint(
                        geo_types::MultiPoint(pts)
                    )
                )
            },
            GeomType::Linestring => {
                Ok(
                    geo_types::Geometry::LineString(
                        geo_line_from_sfcgal(unsafe { self.0.as_ref() })?
                    )
                )
            },
            GeomType::Multilinestring => {
                let ngeoms = unsafe {
                    sfcgal_geometry_collection_num_geometries(self.0.as_ref())
                };
                let mut lines = Vec::with_capacity(ngeoms);
                for i in 0..ngeoms {
                    let geom = unsafe { sfcgal_geometry_collection_geometry_n(self.0.as_ref(), i) };
                    lines.push(geo_line_from_sfcgal(geom)?);
                }
                Ok(
                    geo_types::Geometry::MultiLineString(
                        geo_types::MultiLineString(lines)
                    )
                )
            },
            GeomType::Polygon => {
                let nrings = unsafe { sfcgal_polygon_num_interior_rings(self.0.as_ref()) };
                let exterior_sfcgal = unsafe { sfcgal_polygon_exterior_ring(self.0.as_ref()) };
                let exterior_geo = geo_line_from_sfcgal(exterior_sfcgal)?;
                let mut interiors_geo = Vec::with_capacity(nrings);
                for i in 0..nrings {
                    let line_sfcgal = unsafe {
                        sfcgal_polygon_interior_ring_n(self.0.as_ref(), i)
                    };
                    interiors_geo.push(geo_line_from_sfcgal(line_sfcgal)?);
                }

                Ok(
                    geo_types::Geometry::Polygon(
                        geo_types::Polygon::new(exterior_geo, interiors_geo)
                    )
                )
            }
            GeomType::Multipolygon => {
                let ngeoms = unsafe {
                    sfcgal_geometry_collection_num_geometries(self.0.as_ref())
                };
                let mut vec_polygons = Vec::with_capacity(ngeoms);
                for i in 0..ngeoms {
                    let _polyg = unsafe { sfcgal_geometry_collection_geometry_n(self.0.as_ref(), i) };
                    let nrings = unsafe { sfcgal_polygon_num_interior_rings(_polyg) };
                    let exterior_sfcgal = unsafe { sfcgal_polygon_exterior_ring(_polyg) };
                    let exterior_geo = geo_line_from_sfcgal(exterior_sfcgal)?;
                    let mut interiors_geo = Vec::with_capacity(nrings);
                    for j in 0..nrings {
                        let line_sfcgal = unsafe {
                            sfcgal_polygon_interior_ring_n(_polyg, j)
                        };
                        interiors_geo.push(geo_line_from_sfcgal(line_sfcgal)?);
                    }
                    vec_polygons.push(geo_types::Polygon::new(exterior_geo, interiors_geo));
                }

                Ok(
                    geo_types::Geometry::MultiPolygon(
                        geo_types::MultiPolygon(vec_polygons)
                    )
                )
            }
            _ => unimplemented!()
        }
    }
}

fn geo_line_from_sfcgal(sfcgal_geom: *const sfcgal_geometry_t) -> Result<geo_types::LineString<f64>> {
    let n_points = unsafe { sfcgal_linestring_num_points(sfcgal_geom) };
    let mut v_points = Vec::with_capacity(n_points);
    for i in 0..n_points {
        let pt_sfcgal = unsafe { sfcgal_linestring_point_n(sfcgal_geom, i) };
        check_null_geom(pt_sfcgal)?;
        let pt_geom = geo_point_from_sfcgal(pt_sfcgal);
        v_points.push(pt_geom);
    }
    Ok(geo_types::LineString::from(v_points))
}
fn geo_point_from_sfcgal(geom: *const sfcgal_geometry_t) -> geo_types::Point<f64> {
    let x = unsafe { sfcgal_point_x(geom) };
    let y = unsafe { sfcgal_point_y(geom) };
    geo_types::Point::new(x, y)
}

impl ToSfcgal for geo_types::Point<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let geom = unsafe { sfcgal_point_create_from_xy(self.x(), self.y()) };
        unsafe { SFCGeometry::new_from_raw(geom) }
    }
}

impl ToSfcgal for geo_types::MultiPoint<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_geom = unsafe { sfcgal_multi_point_create() };
        check_null_geom(out_geom)?;
        let &geo_types::MultiPoint(ref point_list) = self;
        for point in point_list.iter() {
            let geom = unsafe {
                sfcgal_point_create_from_xy(point.x(), point.y())
            };
            check_null_geom(geom)?;
            unsafe {
                sfcgal_geometry_collection_add_geometry(out_geom, geom)
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_geom) }
    }
}

impl ToSfcgal for geo_types::Line<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_linestring = unsafe { sfcgal_linestring_create() };
        check_null_geom(out_linestring)?;
        let start = unsafe {
            sfcgal_point_create_from_xy(self.start.x, self.start.y)
        };
        let end = unsafe {
            sfcgal_point_create_from_xy(self.end.x, self.end.y)
        };
        check_null_geom(start)?;
        check_null_geom(end)?;
        unsafe {
            sfcgal_linestring_add_point(out_linestring, start);
            sfcgal_linestring_add_point(out_linestring, end);
            SFCGeometry::new_from_raw(out_linestring)
        }
    }
}

impl ToSfcgal for geo_types::LineString<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        unsafe { SFCGeometry::new_from_raw(linestring_geo_to_sfcgal(self)?) }
    }
}

impl ToSfcgal for geo_types::MultiLineString<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_multilinestring = unsafe { sfcgal_multi_linestring_create() };
        check_null_geom(out_multilinestring)?;
        let &geo_types::MultiLineString(ref linestring_list) = self;
        for _linestring in linestring_list.into_iter() {
            let out_sfcgal_line = linestring_geo_to_sfcgal(_linestring)?;
            unsafe {
                sfcgal_geometry_collection_add_geometry(out_multilinestring, out_sfcgal_line)
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_multilinestring) }
    }
}

impl ToSfcgal for geo_types::Polygon<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let &geo_types::Polygon{ref exterior, ref interiors} = self;

        let out_polygon = unsafe {
            sfcgal_polygon_create_from_exterior_ring(linestring_geo_to_sfcgal(exterior)?)
        };
        check_null_geom(out_polygon)?;

        for ring in interiors {
            unsafe {
                sfcgal_polygon_add_interior_ring(out_polygon, linestring_geo_to_sfcgal(ring)?)
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_polygon) }
    }
}

impl ToSfcgal for geo_types::MultiPolygon<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_multipolygon = unsafe { sfcgal_multi_polygon_create() };
        let &geo_types::MultiPolygon(ref list_polygons) = self;
        for polygon in list_polygons {
            let &geo_types::Polygon{ref exterior, ref interiors} = polygon;
            let out_polygon = unsafe {
                sfcgal_polygon_create_from_exterior_ring(linestring_geo_to_sfcgal(exterior)?)
            };
            check_null_geom(out_polygon)?;

            for ring in interiors {
                unsafe {
                    sfcgal_polygon_add_interior_ring(out_polygon, linestring_geo_to_sfcgal(ring)?)
                };
            }
            unsafe {
                sfcgal_geometry_collection_add_geometry(out_multipolygon, out_polygon)
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_multipolygon) }
    }
}

impl ToSfcgal for geo_types::GeometryCollection<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_geom_collection = unsafe { sfcgal_geometry_collection_create() };
        let &geo_types::GeometryCollection(ref list_geoms) = self;
        for g_geom in list_geoms {
            let sfcgal_geom = g_geom.to_sfcgal()?;
            unsafe {
                sfcgal_geometry_collection_add_geometry(out_geom_collection, sfcgal_geom.0.as_ptr())
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_geom_collection) }
    }
}


impl ToSfcgal for geo_types::Geometry<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        match *self {
            geo_types::Geometry::Point(ref c) => c.to_sfcgal(),
            geo_types::Geometry::Line(ref c) => c.to_sfcgal(),
            geo_types::Geometry::LineString(ref c) => c.to_sfcgal(),
            geo_types::Geometry::Polygon(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiPoint(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiLineString(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiPolygon(ref c) => c.to_sfcgal(),
            geo_types::Geometry::GeometryCollection(ref c) => c.to_sfcgal()
        }
    }
}

fn linestring_geo_to_sfcgal(geom: &geo_types::LineString<f64>) -> Result<*mut sfcgal_geometry_t> {
    let out_linestring = unsafe { sfcgal_linestring_create() };
    check_null_geom(out_linestring)?;
    let &geo_types::LineString(ref point_list) = geom;
    for coords in point_list.iter() {
        let geom = unsafe {
            sfcgal_point_create_from_xy(coords.x, coords.y)
        };
        check_null_geom(geom)?;
        unsafe {
            sfcgal_linestring_add_point(out_linestring, geom)
        };
    }
    Ok(out_linestring)
}


#[cfg(test)]
mod test {
    use geo_types::{Geometry, Point, MultiPoint, LineString, MultiLineString, Polygon, MultiPolygon};
    use crate::{SFCGeometry, ToSfcgal};
    use super::TryInto;

    #[test]
    fn point_geo_to_sfcgal() {
        let pt = Point::new(0.1, 0.9);
        let pt_sfcgal = pt.to_sfcgal().unwrap();
        let pt: Point<f64> = pt_sfcgal.try_into().unwrap().as_point().unwrap();
        assert_eq!(pt.x(), 0.1);
        assert_eq!(pt.y(), 0.9);
    }

    #[test]
    fn point_sfcgal_try_into_geo() {
        let pt_sfcgal = SFCGeometry::new("POINT(0.1 0.9)").unwrap();
        let pt: Point<f64> = pt_sfcgal.try_into().unwrap().as_point().unwrap();
        assert_ulps_eq!(pt.x(), 0.1);
        assert_ulps_eq!(pt.y(), 0.9);
    }

    #[test]
    fn multipoint_geo_to_sfcgal_to_geo() {
        let multipt = MultiPoint::from(vec![
            Point::new(0., 0.),
            Point::new(1., 1.),
        ]);
        let pt_sfcgal = multipt.to_sfcgal().unwrap();
        let mpt: MultiPoint<f64> = pt_sfcgal.try_into().unwrap().as_multipoint().unwrap();
        assert_eq!(mpt.0[0].x(), 0.);
        assert_eq!(mpt.0[0].y(), 0.);
        assert_eq!(mpt.0[1].x(), 1.);
        assert_eq!(mpt.0[1].y(), 1.);
    }

    #[test]
    fn linestring_geo_to_sfcgal_to_geo() {
        let linestring = LineString::from(vec![
            Point::new(0., 0.),
            Point::new(1., 1.),
        ]);
        let line_sfcgal = linestring.to_sfcgal().unwrap();
        let linestring_geo: LineString<f64> = line_sfcgal.try_into().unwrap().as_linestring().unwrap();
        assert_eq!(linestring_geo.0[0].x, 0.);
        assert_eq!(linestring_geo.0[0].y, 0.);
        assert_eq!(linestring_geo.0[1].x, 1.);
        assert_eq!(linestring_geo.0[1].y, 1.);
    }

    #[test]
    fn multilinestring_geo_to_sfcgal_to_geo() {
        let multilinestring = MultiLineString::from(LineString::from(vec![
            Point::new(0., 0.),
            Point::new(1., 1.),
        ]));
        let mls_sfcgal = multilinestring.to_sfcgal().unwrap();
        let mls: MultiLineString<f64> = mls_sfcgal.try_into().unwrap().as_multilinestring().unwrap();
        assert_eq!(mls.0[0].0[0].x, 0.);
        assert_eq!(mls.0[0].0[0].y, 0.);
        assert_eq!(mls.0[0].0[1].x, 1.);
        assert_eq!(mls.0[0].0[1].y, 1.);
    }

    #[test]
    fn polygon_geo_to_sfcgal_to_geo() {
        let polygon = Polygon::new(
            LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)]),
            vec![LineString::from(
                vec![(0.1, 0.1), (0.1, 0.9,), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)])]);
        let poly_sfcgal = polygon.to_sfcgal().unwrap();
        let polyg: Polygon<f64> = poly_sfcgal.try_into().unwrap().as_polygon().unwrap();

        assert_eq!(
            polyg.exterior,
            LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)]));
        assert_eq!(polyg.interiors[0].0[0].x, 0.1);
        assert_eq!(polyg.interiors[0].0[0].y, 0.1);
        assert_eq!(polyg.interiors[0].0[2].x, 0.9);
        assert_eq!(polyg.interiors[0].0[2].y, 0.9);
        assert_eq!(polyg.interiors[0].0[3].x, 0.9);
        assert_eq!(polyg.interiors[0].0[3].y, 0.1);
    }

    #[test]
    fn multipolygon_geo_to_sfcgal_to_geo() {
        let multipolygon = MultiPolygon(vec![
            Polygon::new(
                LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)]),
                vec![LineString::from(
                    vec![(0.1, 0.1), (0.1, 0.9,), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)])]
            ),
        ]);
        let mutlipolygon_sfcgal = multipolygon.to_sfcgal().unwrap();
        let mpg: MultiPolygon<f64> = mutlipolygon_sfcgal.try_into().unwrap().as_multipolygon().unwrap();

        assert_eq!(
            mpg.0[0].exterior,
            LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)]));
        assert_eq!(
            mpg.0[0].interiors[0],
            LineString::from(
                vec![(0.1, 0.1), (0.1, 0.9,), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)]));
    }
}
