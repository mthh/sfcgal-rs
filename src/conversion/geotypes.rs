use crate::conversion::coords::{CoordSeq, CoordType, ToSFCGALGeom};
use crate::{
    utils::check_null_geom, GeomType, Point2d, Result, SFCGeometry, ToCoordinates, ToSFCGAL,
};
use anyhow::Error;
use sfcgal_sys::{
    sfcgal_geometry_collection_add_geometry, sfcgal_geometry_collection_create,
    sfcgal_geometry_collection_geometry_n, sfcgal_geometry_collection_num_geometries,
    sfcgal_geometry_t, sfcgal_linestring_add_point, sfcgal_linestring_create,
    sfcgal_linestring_num_points, sfcgal_linestring_point_n, sfcgal_multi_linestring_create,
    sfcgal_multi_point_create, sfcgal_multi_polygon_create, sfcgal_point_create_from_xy,
    sfcgal_point_x, sfcgal_point_y, sfcgal_polygon_add_interior_ring,
    sfcgal_polygon_create_from_exterior_ring, sfcgal_polygon_exterior_ring,
    sfcgal_polygon_interior_ring_n, sfcgal_polygon_num_interior_rings, sfcgal_triangle_create,
    sfcgal_triangle_set_vertex_from_xy,
};
use std::convert::Into;
use std::iter::FromIterator;

/// Conversion from [`SFCGeometry`] (implemented on [geo-types](https://docs.rs/geo-types/) geometries)
///
/// [`SFCGeometry`]: struct.SFCGeometry.html
pub trait TryInto<T> {
    type Err;
    fn try_into(self) -> Result<T>;
}

impl CoordType for geo_types::Point<f64> {}
impl CoordType for geo_types::Coordinate<f64> {}

impl ToSFCGALGeom for geo_types::Point<f64> {
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t> {
        let g = unsafe { sfcgal_point_create_from_xy(self.x(), self.y()) };
        check_null_geom(g)?;
        Ok(g)
    }
}

impl ToSFCGALGeom for geo_types::Coordinate<f64> {
    fn to_sfcgeometry(&self) -> Result<*mut sfcgal_geometry_t> {
        let g = unsafe { sfcgal_point_create_from_xy(self.x, self.y) };
        check_null_geom(g)?;
        Ok(g)
    }
}

/// Implements conversion from CoordSeq to geo_types::Geometry
/// (better use TryInto<geo_types::Geometry> for SFCGeometry if the intend
/// is to convert SFCGAL Geometries to geo_types ones)
impl TryInto<geo_types::Geometry<f64>> for CoordSeq<Point2d> {
    type Err = Error;

    fn try_into(self) -> Result<geo_types::Geometry<f64>> {
        match self {
                CoordSeq::Point(pt) => Ok(geo_types::Point(pt.into()).into()),
                CoordSeq::Multipoint(pts) => Ok(geo_types::MultiPoint::from_iter(pts.into_iter()).into()),
                CoordSeq::Linestring(pts) => Ok(geo_types::LineString::from_iter(pts.into_iter()).into()),
                CoordSeq::Multilinestring(lines) => {
                    Ok(geo_types::MultiLineString(
                        lines.into_iter()
                            .map(geo_types::LineString::from)
                            .collect()
                    ).into())
                },
                CoordSeq::Polygon(rings) => {
                    let mut it = rings.into_iter();
                    let exterior = geo_types::LineString::from(it.next().unwrap());
                    let interiors = it.map(geo_types::LineString::from).collect::<Vec<geo_types::LineString<f64>>>();
                    Ok(geo_types::Polygon::new(exterior, interiors).into())
                },
                CoordSeq::Multipolygon(polygons) => {
                    let polys = polygons.into_iter().map(|p| {
                        let a: geo_types::Geometry<f64> = CoordSeq::Polygon(p).try_into()?;
                        if let Some(poly) = a.into_polygon() {
                            Ok(poly)
                        } else {
                            Err(format_err!("Error while building geo_types::MultiPolygon"))
                        }
                    }).collect::<Result<Vec<geo_types::Polygon<f64>>>>()?;
                    Ok(geo_types::MultiPolygon(polys).into())
                },
                CoordSeq::Geometrycollection(collection) => {
                    Ok(geo_types::Geometry::GeometryCollection(
                        geo_types::GeometryCollection(collection
                            .into_iter()
                            .map(|g| g.try_into())
                            .collect::<Result<Vec<geo_types::Geometry<f64>>>>()?
                        )
                    ))
                },
                _ => Err(
                    format_err!(
                        "Conversion from CoordSeq variants `Solid`, `Multisolid`, `Triangulatedsurface` and `Polyhedralsurface` are not yet implemented!"))
            }
    }
}

/// Implements faillible conversion from SFCGeometry to geo_types::Geometry.
///
/// This is notably faillible because some types of [`SFCGeometry`] like GeoTypes::Polyhedralsurface
/// don't have equivalents in geo_types::Geometry.
/// Please note that geo_types Coordinate and Point primitives are 2d only, so
/// every information about z coordinate (if any) won't be taken into account.
impl TryInto<geo_types::Geometry<f64>> for SFCGeometry {
    type Err = Error;

    fn try_into(self) -> Result<geo_types::Geometry<f64>> {
        match self._type()? {
            GeomType::Point => {
                let c = self.to_coordinates::<Point2d>()?;
                let p: geo_types::Point<f64> = match c {
                    CoordSeq::Point(pt) => pt.into(),
                    _ => unimplemented!(),
                };
                Ok(geo_types::Geometry::Point(p))
            }
            GeomType::Multipoint => {
                let c = self.to_coordinates::<Point2d>()?;
                let p: geo_types::MultiPoint<f64> = match c {
                    CoordSeq::Multipoint(pts) => pts.into(),
                    _ => unimplemented!(),
                };
                Ok(geo_types::Geometry::MultiPoint(p))
            }
            GeomType::Linestring => Ok(geo_types::Geometry::LineString(geo_line_from_sfcgal(
                unsafe { self.c_geom.as_ref() },
            )?)),
            GeomType::Multilinestring => {
                let ngeoms =
                    unsafe { sfcgal_geometry_collection_num_geometries(self.c_geom.as_ref()) };
                let mut lines = Vec::with_capacity(ngeoms);
                for i in 0..ngeoms {
                    let geom =
                        unsafe { sfcgal_geometry_collection_geometry_n(self.c_geom.as_ref(), i) };
                    lines.push(geo_line_from_sfcgal(geom)?);
                }
                Ok(geo_types::Geometry::MultiLineString(
                    geo_types::MultiLineString(lines),
                ))
            }
            GeomType::Polygon => {
                let nrings = unsafe { sfcgal_polygon_num_interior_rings(self.c_geom.as_ref()) };
                let exterior_sfcgal = unsafe { sfcgal_polygon_exterior_ring(self.c_geom.as_ref()) };
                let exterior_geo = geo_line_from_sfcgal(exterior_sfcgal)?;
                let mut interiors_geo = Vec::with_capacity(nrings);
                for i in 0..nrings {
                    let line_sfcgal =
                        unsafe { sfcgal_polygon_interior_ring_n(self.c_geom.as_ref(), i) };
                    interiors_geo.push(geo_line_from_sfcgal(line_sfcgal)?);
                }

                Ok(geo_types::Geometry::Polygon(geo_types::Polygon::new(
                    exterior_geo,
                    interiors_geo,
                )))
            }
            GeomType::Multipolygon => {
                let ngeoms =
                    unsafe { sfcgal_geometry_collection_num_geometries(self.c_geom.as_ref()) };
                let mut vec_polygons = Vec::with_capacity(ngeoms);
                for i in 0..ngeoms {
                    let _polyg =
                        unsafe { sfcgal_geometry_collection_geometry_n(self.c_geom.as_ref(), i) };
                    let nrings = unsafe { sfcgal_polygon_num_interior_rings(_polyg) };
                    let exterior_sfcgal = unsafe { sfcgal_polygon_exterior_ring(_polyg) };
                    let exterior_geo = geo_line_from_sfcgal(exterior_sfcgal)?;
                    let mut interiors_geo = Vec::with_capacity(nrings);
                    for j in 0..nrings {
                        let line_sfcgal = unsafe { sfcgal_polygon_interior_ring_n(_polyg, j) };
                        interiors_geo.push(geo_line_from_sfcgal(line_sfcgal)?);
                    }
                    vec_polygons.push(geo_types::Polygon::new(exterior_geo, interiors_geo));
                }

                Ok(geo_types::MultiPolygon(vec_polygons).into())
            }
            GeomType::Geometrycollection => {
                let c = self.to_coordinates::<Point2d>()?;
                let p = match c {
                    CoordSeq::Geometrycollection(g) => g
                        .into_iter()
                        .map(|g| g.try_into())
                        .collect::<Result<Vec<geo_types::Geometry<f64>>>>()?,
                    _ => unimplemented!(),
                };
                Ok(geo_types::Geometry::GeometryCollection(
                    geo_types::GeometryCollection(p),
                ))
            }
            _ => Err(format_err!(
                "Conversion from SFCGeometry of type `Triangle`, `Solid`, `Multisolid`, \
                 `Triangulatedsurface` and `Polyhedralsurface` \
                 to geo_types::Geometry are not yet implemented!"
            )),
        }
    }
}

fn geo_line_from_sfcgal(
    sfcgal_geom: *const sfcgal_geometry_t,
) -> Result<geo_types::LineString<f64>> {
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

/// Create a `SFCGeometry` from a geo-types Point
impl ToSFCGAL for geo_types::Point<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        unsafe { SFCGeometry::new_from_raw(sfcgal_point_create_from_xy(self.x(), self.y()), true) }
    }
}

/// Create a `SFCGeometry` from a geo-types MultiPoint
impl ToSFCGAL for geo_types::MultiPoint<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        make_sfcgal_multi_geom!(
            sfcgal_multi_point_create(),
            self.0
                .iter()
                .map(|pt| pt.to_sfcgeometry())
                .collect::<Result<Vec<_>>>()?
        )
    }
}

/// Create a `SFCGeometry` from a geo-types Line
impl ToSFCGAL for geo_types::Line<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_linestring = unsafe { sfcgal_linestring_create() };
        check_null_geom(out_linestring)?;
        let start = unsafe { sfcgal_point_create_from_xy(self.start.x, self.start.y) };
        let end = unsafe { sfcgal_point_create_from_xy(self.end.x, self.end.y) };
        check_null_geom(start)?;
        check_null_geom(end)?;
        unsafe {
            sfcgal_linestring_add_point(out_linestring, start);
            sfcgal_linestring_add_point(out_linestring, end);
            SFCGeometry::new_from_raw(out_linestring, true)
        }
    }
}
/// Create a `SFCGeometry` from a geo-types LineString
impl ToSFCGAL for geo_types::LineString<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let line = (&self.0).to_sfcgeometry()?;
        unsafe { SFCGeometry::new_from_raw(line, true) }
    }
}

/// Create a `SFCGeometry` from a geo-types MultiLineString
impl ToSFCGAL for geo_types::MultiLineString<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        make_sfcgal_multi_geom!(
            sfcgal_multi_linestring_create(),
            self.0
                .iter()
                .map(|line| line.0.to_sfcgeometry())
                .collect::<Result<Vec<_>>>()?
        )
    }
}

/// Create a `SFCGeometry` from a geo-types Triangle
impl ToSFCGAL for geo_types::Triangle<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        let out_triangle = unsafe { sfcgal_triangle_create() };
        check_null_geom(out_triangle)?;
        let &geo_types::Triangle(ref c0, ref c1, ref c2) = self;
        unsafe {
            sfcgal_triangle_set_vertex_from_xy(out_triangle, 0, c0.x, c0.y);
            sfcgal_triangle_set_vertex_from_xy(out_triangle, 1, c1.x, c1.y);
            sfcgal_triangle_set_vertex_from_xy(out_triangle, 2, c2.x, c2.y);
            SFCGeometry::new_from_raw(out_triangle, true)
        }
    }
}

fn geo_polygon_to_sfcgal<T>(
    exterior: &Vec<T>,
    interiors: &[geo_types::LineString<f64>],
) -> Result<*mut sfcgal_geometry_t>
where
    T: ToSFCGALGeom + CoordType,
{
    let out_polygon =
        unsafe { sfcgal_polygon_create_from_exterior_ring(exterior.to_sfcgeometry()?) };
    check_null_geom(out_polygon)?;
    for ring in interiors.iter() {
        unsafe { sfcgal_polygon_add_interior_ring(out_polygon, ring.0.to_sfcgeometry()?) };
    }
    Ok(out_polygon)
}

/// Create a `SFCGeometry` from a geo-types Polygon
impl ToSFCGAL for geo_types::Polygon<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        // let geo_types::Polygon{exterior, interiors} = self;
        let (exterior, interiors) = (self.exterior(), self.interiors());
        let out_polygon = geo_polygon_to_sfcgal(&exterior.0, &interiors)?;
        unsafe { SFCGeometry::new_from_raw(out_polygon, true) }
    }
}

/// Create a `SFCGeometry` from a geo-types MultiPolygon
impl ToSFCGAL for geo_types::MultiPolygon<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        make_sfcgal_multi_geom!(
            sfcgal_multi_polygon_create(),
            self.0
                .iter()
                .map(|polygon| {
                    // let geo_types::Polygon{ref exterior, ref interiors} = polygon;
                    let (exterior, interiors) = (polygon.exterior(), polygon.interiors());
                    geo_polygon_to_sfcgal(&exterior.0, &interiors)
                })
                .collect::<Result<Vec<_>>>()?
        )
    }
}

/// Create a `SFCGeometry` from a geo-types GeometryCollection
impl ToSFCGAL for geo_types::GeometryCollection<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        make_sfcgal_multi_geom!(
            sfcgal_geometry_collection_create(),
            self.0
                .iter()
                .map(|geom| {
                    let mut _g = geom.to_sfcgal()?;
                    _g.owned = false;
                    Ok(_g.c_geom.as_ptr())
                })
                .collect::<Result<Vec<_>>>()?
        )
    }
}

/// Create a `SFCGeometry` from any geo-type Geometry
impl ToSFCGAL for geo_types::Geometry<f64> {
    fn to_sfcgal(&self) -> Result<SFCGeometry> {
        match *self {
            geo_types::Geometry::Point(ref c) => c.to_sfcgal(),
            geo_types::Geometry::Line(ref c) => c.to_sfcgal(),
            geo_types::Geometry::LineString(ref c) => c.to_sfcgal(),
            geo_types::Geometry::Polygon(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiPoint(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiLineString(ref c) => c.to_sfcgal(),
            geo_types::Geometry::MultiPolygon(ref c) => c.to_sfcgal(),
            geo_types::Geometry::GeometryCollection(ref c) => c.to_sfcgal(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TryInto;
    use crate::{GeomType, SFCGeometry, ToSFCGAL};
    use geo_types::{
        Coordinate, LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, Triangle,
    };

    #[test]
    fn point_geo_to_sfcgal_to_geo() {
        let pt = Point::new(0.1, 0.9);
        let pt_sfcgal = pt.to_sfcgal().unwrap();
        assert!(pt_sfcgal.is_valid().unwrap());
        let pt: Point<f64> = pt_sfcgal.try_into().unwrap().into_point().unwrap();
        assert_eq!(pt.x(), 0.1);
        assert_eq!(pt.y(), 0.9);
    }

    #[test]
    fn point_sfcgal_try_into_geo() {
        let pt_sfcgal = SFCGeometry::new("POINT(0.1 0.9)").unwrap();
        let pt: Point<f64> = pt_sfcgal.try_into().unwrap().into_point().unwrap();
        assert_ulps_eq!(pt.x(), 0.1);
        assert_ulps_eq!(pt.y(), 0.9);
    }

    #[test]
    fn multipoint_geo_to_sfcgal_to_geo() {
        let multipt = MultiPoint::from(vec![Point::new(0., 0.), Point::new(1., 1.)]);
        let mpt_sfcgal = multipt.to_sfcgal().unwrap();
        assert!(mpt_sfcgal.is_valid().unwrap());
        let mpt: MultiPoint<f64> = mpt_sfcgal.try_into().unwrap().into_multi_point().unwrap();
        assert_eq!(mpt.0[0].x(), 0.);
        assert_eq!(mpt.0[0].y(), 0.);
        assert_eq!(mpt.0[1].x(), 1.);
        assert_eq!(mpt.0[1].y(), 1.);
    }

    #[test]
    fn linestring_geo_to_sfcgal_to_geo() {
        let linestring = LineString::from(vec![Point::new(0., 0.), Point::new(1., 1.)]);
        let line_sfcgal = linestring.to_sfcgal().unwrap();
        assert!(line_sfcgal.is_valid().unwrap());
        let linestring_geo: LineString<f64> =
            line_sfcgal.try_into().unwrap().into_line_string().unwrap();
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
        assert!(mls_sfcgal.is_valid().unwrap());
        let mls: MultiLineString<f64> = mls_sfcgal
            .try_into()
            .unwrap()
            .into_multi_line_string()
            .unwrap();
        assert_eq!(mls.0[0].0[0].x, 0.);
        assert_eq!(mls.0[0].0[0].y, 0.);
        assert_eq!(mls.0[0].0[1].x, 1.);
        assert_eq!(mls.0[0].0[1].y, 1.);
    }

    #[test]
    fn triangle_geo_to_sfcgal_to_geo() {
        let tri = Triangle(
            Coordinate::from((0., 0.)),
            Coordinate::from((1., 0.)),
            Coordinate::from((0.5, 1.)),
        );
        let tri_sfcgal = tri.to_sfcgal().unwrap();
        assert!(tri_sfcgal.is_valid().unwrap());
        assert_eq!(tri_sfcgal._type().unwrap(), GeomType::Triangle);
        let coords: Result<geo_types::Geometry<f64>, _> = tri_sfcgal.try_into();
        assert_eq!(
            coords.err().unwrap().to_string(),
            "Conversion from SFCGeometry of type `Triangle`, `Solid`, `Multisolid`, \
            `Triangulatedsurface` and `Polyhedralsurface` to geo_types::Geometry are not yet implemented!",
        )
    }

    #[test]
    fn polygon_geo_to_sfcgal_to_geo() {
        let polygon = Polygon::new(
            LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.), (0., 0.)]),
            vec![LineString::from(vec![
                (0.1, 0.1),
                (0.1, 0.9),
                (0.9, 0.9),
                (0.9, 0.1),
                (0.1, 0.1),
            ])],
        );
        let poly_sfcgal = polygon.to_sfcgal().unwrap();
        let polyg: Polygon<f64> = poly_sfcgal.try_into().unwrap().into_polygon().unwrap();
        let interiors = polyg.interiors();
        assert_eq!(
            polyg.exterior(),
            &LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)])
        );
        assert_eq!(interiors[0].0[0].x, 0.1);
        assert_eq!(interiors[0].0[0].y, 0.1);
        assert_eq!(interiors[0].0[2].x, 0.9);
        assert_eq!(interiors[0].0[2].y, 0.9);
        assert_eq!(interiors[0].0[3].x, 0.9);
        assert_eq!(interiors[0].0[3].y, 0.1);
    }

    #[test]
    fn multipolygon_geo_to_sfcgal_to_geo() {
        let multipolygon = MultiPolygon(vec![Polygon::new(
            LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.), (0., 0.)]),
            vec![LineString::from(vec![
                (0.1, 0.1),
                (0.1, 0.9),
                (0.9, 0.9),
                (0.9, 0.1),
                (0.1, 0.1),
            ])],
        )]);
        let mutlipolygon_sfcgal = multipolygon.to_sfcgal().unwrap();
        let mpg: MultiPolygon<f64> = mutlipolygon_sfcgal
            .try_into()
            .unwrap()
            .into_multi_polygon()
            .unwrap();

        assert_eq!(
            mpg.0[0].exterior(),
            &LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)])
        );
        assert_eq!(
            mpg.0[0].interiors()[0],
            LineString::from(vec![
                (0.1, 0.1),
                (0.1, 0.9,),
                (0.9, 0.9),
                (0.9, 0.1),
                (0.1, 0.1)
            ])
        );
    }

    #[test]
    fn geometrycollection_sfcgal_to_geo_to_sfcgal() {
        let input_wkt = "GEOMETRYCOLLECTION(POINT(4.0 6.0),LINESTRING(4.0 6.0,7.0 10.0))";
        let gc_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let gc: geo_types::Geometry<f64> = gc_sfcgal.try_into().unwrap();
        if let geo_types::Geometry::GeometryCollection(_gc) = &gc {
            assert_eq!(Point::new(4., 6.), _gc.0[0].clone().into_point().unwrap(),);
            assert_eq!(
                LineString::from(vec![(4., 6.), (7., 10.)]),
                _gc.0[1].clone().into_line_string().unwrap(),
            );
            let gc_sfcgal = _gc.to_sfcgal().unwrap();
            assert_eq!(input_wkt, gc_sfcgal.to_wkt_decim(1).unwrap());
        } else {
            panic!("Error while deconstructing geometrycollection");
        }
    }
}
