use geojson::{Value as GeometryValue};
use crate::{ToSFCGAL, CoordSeq, Result, SFCGeometry, ToCoordinates, Point2d, Point3d};
use crate::conversion::coords::{ CoordType, ToSFCGALGeom, FromSFCGALGeom};
use failure::Error;


/// Conversion from [`SFCGeometry`] (implemented on [geo-types](https://docs.rs/geo-types/) geometries)
///
/// [`SFCGeometry`]: struct.SFCGeometry.html
pub trait TryIntoCoords<T> {
    type Err;
    fn try_into(&self) -> Result<T>;
}

/// Convert coordinates to `sfcgal_geometry_t`.
pub trait FromSlice {
    fn from_slice(pt: &[f64]) -> Result<Self> where Self: std::marker::Sized;
}

/// Convert coordinates to `sfcgal_geometry_t`.
pub trait ToVec {
    fn to_vec(&self) -> Vec<f64>;
}

impl ToVec for Point2d {
    fn to_vec(&self) -> Vec<f64> {
        vec![self.0, self.1]
    }
}

impl ToVec for Point3d {
    fn to_vec(&self) -> Vec<f64> {
        vec![self.0, self.1, self.2]
    }
}

impl FromSlice for Point2d {
    fn from_slice(pt: &[f64]) -> Result<Point2d> {
        let mut it = pt.iter();
        Ok((*it.next().unwrap(), *it.next().unwrap()))
    }
}

impl FromSlice for Point3d {
    fn from_slice(pt: &[f64]) -> Result<Point3d> {
        let mut it = pt.iter();
        Ok((*it.next().unwrap(), *it.next().unwrap(), *it.next().unwrap_or_else(|| &(0.0f64))))    }
}

/// Implements conversion from CoordSeq to geo_types::Geometry
/// (better use TryInto<geo_types::Geometry> for SFCGeometry if the intend
/// is to convert SFCGAL Geometries to geo_types ones)
impl<T> TryIntoCoords<CoordSeq<T>> for GeometryValue where T: FromSlice + CoordType {
    type Err = Error;
    fn try_into(&self) -> Result<CoordSeq<T>> where T: FromSlice + CoordType {
        match *self {
            GeometryValue::Point(ref pt) => {
                Ok(CoordSeq::Point(T::from_slice(&pt)?))
            },
            GeometryValue::MultiPoint(ref pts) => {
                let _pts = pts.iter().map(|p| T::from_slice(p)).collect::<Result<Vec<T>>>();
                Ok(CoordSeq::Multipoint(_pts?))
            },
            GeometryValue::LineString(ref pts) => {
                let _pts = pts.iter().map(|p| T::from_slice(p)).collect::<Result<Vec<T>>>();
                Ok(CoordSeq::Linestring(_pts?))
            },
            GeometryValue::MultiLineString(ref lines) => {
                let _pts = lines.iter().map(|pts| {
                    pts.iter().map(|p| T::from_slice(p)).collect()
                }).collect::<Result<Vec<Vec<T>>>>();
                Ok(CoordSeq::Multilinestring(_pts?))
            },
            GeometryValue::Polygon(ref rings) => {
                let _pts = rings.iter().map(|pts| {
                    pts.iter().map(|p| T::from_slice(p)).collect()
                }).collect::<Result<Vec<Vec<T>>>>();
                Ok(CoordSeq::Polygon(_pts?))
            },
            GeometryValue::MultiPolygon(ref polygons) => {
                let _pts = polygons.iter().map(|rings|{
                    rings.iter().map(|pts| {
                        pts.iter().map(|p| T::from_slice(p)).collect()
                    }).collect()
                }).collect::<Result<Vec<Vec<Vec<T>>>>>();
                Ok(CoordSeq::Multipolygon(_pts?))
            },
            GeometryValue::GeometryCollection(ref geoms) => {
                let _geoms = geoms.iter().map(|ref geom| {
                    geom.value.try_into()
                }).collect::<Result<Vec<CoordSeq<T>>>>();
                Ok(CoordSeq::Geometrycollection(_geoms?))
            }
        }
    }
}

/// Conversion from GeoJson to SFCGAL Geometries.
pub trait FromGeoJSON {
    type Err;
    fn from_geojson<T: FromSlice + CoordType + ToSFCGALGeom>(geom: &GeometryValue) -> Result<SFCGeometry>;
}

/// Conversion from SFCGAL Geometries to GeoJson.
pub trait ToGeoJSON {
    type Err;
    fn to_geojson<T: FromSlice + CoordType + ToVec + FromSFCGALGeom>(&self) -> Result<GeometryValue>;
}

/// Conversion from SFCGAL Geometries to GeoJson.
///
/// Allows to choose if coordinates of constructed geojson have to be 2d or 3d.
/// ``` rust
/// use sfcgal::{SFCGeometry, ToGeoJSON};
/// type Point3d = (f64, f64, f64);
///
/// let input_wkt = "POINT(0.1 0.9 1.0)";
/// let pt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
/// let pt_geojson = pt_sfcgal.to_geojson::<Point3d>().unwrap();
/// ```
impl ToGeoJSON for SFCGeometry {
    type Err = Error;
    fn to_geojson<T: FromSlice + CoordType + ToVec + FromSFCGALGeom>(&self) -> Result<GeometryValue> {
        let cs: CoordSeq<T> = self.to_coordinates()?;
        cs.to_geojson::<T>()
    }
}

impl<U> ToGeoJSON for CoordSeq<U> where U: FromSlice + CoordType + ToVec + FromSFCGALGeom {
    type Err = Error;
    fn to_geojson<T: FromSlice + CoordType + ToVec + FromSFCGALGeom>(&self) -> Result<GeometryValue> {
        match self {
            CoordSeq::Point(pt) => Ok(GeometryValue::Point(pt.to_vec())),
            CoordSeq::Multipoint(pts) => Ok(GeometryValue::MultiPoint(pts.iter().map(|p| p.to_vec()).collect())),
            CoordSeq::Linestring(pts) => Ok(GeometryValue::LineString(pts.iter().map(|p| p.to_vec()).collect())),
            CoordSeq::Multilinestring(lines) => Ok(GeometryValue::MultiLineString(lines.iter().map(|pts| pts.iter().map(|p| p.to_vec()).collect()).collect())),
            CoordSeq::Polygon(rings) => Ok(GeometryValue::Polygon(rings.iter().map(|pts| pts.iter().map(|p| p.to_vec()).collect()).collect())),
            CoordSeq::Multipolygon(polygons) => Ok(
                GeometryValue::MultiPolygon(polygons.iter().map(|rings| rings.iter().map(|pts| pts.iter().map(|p| p.to_vec()).collect()).collect()).collect())),
            CoordSeq::Geometrycollection(geoms) => Ok(
                GeometryValue::GeometryCollection(
                    geoms.iter()
                        .map(|geom|{
                            Ok(geojson::Geometry {
                                bbox: None,
                                value:  geom.to_geojson::<T>()?,
                                foreign_members: None,
                            })
                        }).collect::<Result<Vec<_>>>()?
                    )
                ),
            _ => unimplemented!(),
        }
    }
}

/// Conversion from GeoJson to SFCGAL geometries.
///
/// Allows to choose if the constructed SFCGAL geometries will be 2d or 3d.
/// ``` rust
/// use sfcgal::{SFCGeometry, FromGeoJSON};
/// type Point2d = (f64, f64);
///
/// let geom = geojson::Geometry {
///     bbox: None,
///     value: geojson::Value::LineString(vec![vec![101.0, 0.0], vec![102.0, 1.0]]),
///     foreign_members: None,
/// };
/// let line_sfcgal = SFCGeometry::from_geojson::<Point2d>(&geom.value).unwrap();
/// ```
impl FromGeoJSON for SFCGeometry {
    type Err = Error;
    fn from_geojson<T: FromSlice + CoordType + ToSFCGALGeom>(geom: &GeometryValue) -> Result<SFCGeometry> {
        let coords: CoordSeq<T> = geom.try_into()?;
        coords.to_sfcgal()
    }
}


#[cfg(test)]
mod tests {
    use crate::*;
    use super::{Point2d, Point3d, ToGeoJSON};

    #[test]
    fn point_2d_sfcgal_to_geojson_to_sfcgal() {
        let input_wkt = "POINT(0.1 0.9)";
        let pt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let pt_geojson = pt_sfcgal.to_geojson::<Point2d>().unwrap();
        let pt_sfcgal_new = SFCGeometry::from_geojson::<Point2d>(&pt_geojson).unwrap();
        assert_eq!(pt_sfcgal.to_wkt_decim(1).unwrap(), pt_sfcgal_new.to_wkt_decim(1).unwrap());
    }

    #[test]
    fn point_3d_sfcgal_to_geojson_to_sfcgal() {
        let input_wkt = "POINT(0.1 0.9 1.0)";
        let pt_sfcgal = SFCGeometry::new(input_wkt).unwrap();
        let pt_geojson = pt_sfcgal.to_geojson::<Point3d>().unwrap();
        let pt_sfcgal_new = SFCGeometry::from_geojson::<Point3d>(&pt_geojson).unwrap();
        assert_eq!(pt_sfcgal.to_wkt_decim(1).unwrap(), pt_sfcgal_new.to_wkt_decim(1).unwrap());
    }
}
