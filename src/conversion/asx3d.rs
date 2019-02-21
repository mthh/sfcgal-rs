use crate::{CoordSeq, Result, Point3d};

/// Contains the formatted node element of a geometry in X3D xml
/// and its min / max coordinates in each 3 dimensions.
#[derive(Debug)]
pub struct X3dString {
    pub geometry: String,
    pub bbox: Bbox
}


/// Representation as [`X3dString`] (implemented on some [`CoordSeq`] variants
/// allowing to use it on some [`SFCGeometry`]).
///
/// [`X3dString`]: struct.X3dString.html
/// [`CoordSeq`]: struct.CoordSeq.html
/// [`SFCGeometry`]: struct.SFCGeometry.html
pub trait AsX3d {
    fn as_x3d<S>(self, precision: usize, flip: bool, id: S) -> Result<X3dString> where S: Into<Option<String>>;
}

/// Container for the minimum and maximum values in x, y and z dimensions.
#[derive(Debug)]
pub struct Bbox {
    xmin: f64,
    ymin: f64,
    zmin: f64,
    xmax: f64,
    ymax: f64,
    zmax: f64,
}

impl Bbox {
    pub fn new() -> Self {
        Bbox {
            xmin: std::f64::INFINITY,
            ymin: std::f64::INFINITY,
            zmin: std::f64::INFINITY,
            xmax: std::f64::NEG_INFINITY,
            ymax: std::f64::NEG_INFINITY,
            zmax: std::f64::NEG_INFINITY,
        }
    }

    pub fn process_point(&mut self, pt: &(f64, f64, f64)) {
        if pt.0 < self.xmin {
            self.xmin = pt.0;
        } else if pt.0 > self.xmax {
            self.xmax = pt.0;
        }
        if pt.1 < self.ymin {
            self.ymin = pt.1;
        } else if pt.1 > self.ymax {
            self.ymax = pt.1;
        }
        if pt.2 < self.zmin {
            self.zmin = pt.2;
        } else if pt.2 > self.zmax {
            self.zmax = pt.2;
        }
    }
}

macro_rules! point_string {
    ($geom: expr, $precision: expr, $flip: expr, $bbox: ident) => ({
        $bbox.process_point($geom);
        if $flip == true {
            format!("{:.*} {:.*} {:.*}", $precision, $geom.1, $precision, $geom.0, $precision, $geom.2)
        } else {
            format!("{:.*} {:.*} {:.*}", $precision, $geom.0, $precision, $geom.1, $precision, $geom.2)
        }
    });
}

fn array_point_string(points: &[Point3d], precision: usize, flip: bool, skip_last: bool, bbox: &mut Bbox) -> String {
    points.iter()
        .enumerate()
        .filter_map(|(i, pt)| {
            if skip_last && i == points.len() -1 {
                None
            } else {
                Some(point_string!(pt, precision, flip, bbox))
            }
        })
        .collect::<Vec<String>>()
        .join(" ")
}

fn poly_string(rings: &[Vec<Point3d>], precision: usize, flip: bool, bbox: &mut Bbox) -> String {
    rings.iter()
        .map(|ring|{
            array_point_string(ring, precision, flip, true, bbox)
        })
        .collect::<Vec<String>>()
        .join(" ")
}

/// Representation as [`X3dString`] on [`CoordSeq`] coordinates
/// (as returned by `SFCGeometry.to_coordinates()`). Only implemented for 3d coordinates
/// on : Linestring, Triangulatedsurface and Polyhedralsurface.
///
/// [`X3dString`]: struct.X3dString.html
/// [`CoordSeq`]: struct.CoordSeq.html
impl AsX3d for CoordSeq<Point3d> {
    fn as_x3d<S>(self, precision: usize, flip: bool, id: S) -> Result<X3dString>
        where S: Into<Option<String>>
    {
        let id_obj = match id.into() {
            None => String::default(),
            Some(_id) => format!(" id='{}'", _id)
        };
        let mut b = Bbox::new();
        match self {
            CoordSeq::Linestring(ref pts) => {
                let n_pts = pts.len();
                let s = format!(
                    "<LineSet{} vertexCount='{}'><Coordinate point='{}' /></LineSet>",
                    id_obj,
                    n_pts,
                    array_point_string(&pts, precision, flip, false, &mut b),
                );
                Ok(X3dString {
                    geometry: s,
                    bbox: b,
                })

            },
            CoordSeq::Triangulatedsurface(ref triangles) => {
                let n_triangles = triangles.len();
                let mut s = String::from(format!("<IndexedTriangleSet{} index='", id_obj));
                let mut k = 0;
                for i in 0..n_triangles {
                    s.push_str(&format!("{} {} {}", k, k+1, k+2));
                    if i < n_triangles -1 {
                        s.push_str(" ");
                    }
                    k += 3;
                }
                s.push_str("'><Coordinate point='");
                for i in 0..n_triangles {
                    // No need to skip the last point as closed triangle are only defined
                    // by their 3 point in SFCGAL Geometry / in CoordSeq enum
                    s.push_str(&array_point_string(&triangles[i], precision, flip, false, &mut b));
                    if i < n_triangles - 1 {
                        s.push_str(" ");
                    }
                }
                s.push_str("'/></IndexedTriangleSet>");
                Ok(X3dString {
                    geometry: s,
                    bbox: b,
                })
            },
            CoordSeq::Polyhedralsurface(ref polygons) => {
                let n_geoms = polygons.len();
                let mut s = String::from(format!("<IndexedFaceSet{} coordIndex='", id_obj));
                let mut j = 0;
                polygons
                    .iter()
                    .enumerate()
                    .map(|(i, rings_patch)|{
                        let np: usize = rings_patch.iter().map(|g| g.len()).sum::<usize>() - 1usize;
                        for k in 0..np {
                            if k > 0 { s.push_str(" "); }
                            s.push_str(&format!("{}", (j+ k)));
                        }
                        if i < n_geoms -1 {
                            s.push_str(" -1 ");
                        }
                        j += np;
                    })
                    .for_each(drop);
                s.push_str("'><Coordinate point='");
                s.push_str(&(&polygons.iter()
                    .map(|rings_patch|{
                        poly_string(rings_patch.as_slice(), precision, flip, &mut b)
                    })
                    .collect::<Vec<String>>()
                    .join(" "))
                );
                s.push_str("' /></IndexedFaceSet>");
                Ok(X3dString {
                    geometry: s,
                    bbox: b,
                })
            },
            _ => Err(format_err!("Not implemented for type {:?}", self))
        }
    }
}



#[cfg(test)]
mod tests {
    use crate::{SFCGeometry, ToCoordinates};
    use super::*;

    #[test]
    fn tin_representation() {
        let triangle = SFCGeometry::new("TIN (((0 0 0,0 0 1,0 1 0,0 0 0)), ((0 0 0,0 1 0,1 1 0,0 0 0)))")
            .unwrap()
            .to_coordinates()
            .unwrap();
        let res = triangle.as_x3d(0, false, "tin1".to_string());
        assert_eq!(
            res.unwrap().geometry,
            "<IndexedTriangleSet id='tin1' index='0 1 2 3 4 5'><Coordinate point='0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1 0'/></IndexedTriangleSet>".to_string());
    }

    #[test]
    fn polyhedral_surface_representation() {
        let surface = SFCGeometry::new("POLYHEDRALSURFACE(\
            ((0 0 0, 0 0 1, 0 1 1, 0 1 0, 0 0 0)),\
            ((0 0 0, 0 1 0, 1 1 0, 1 0 0, 0 0 0)),\
            ((0 0 0, 1 0 0, 1 0 1, 0 0 1, 0 0 0)),\
            ((1 1 0, 1 1 1, 1 0 1, 1 0 0, 1 1 0)),\
            ((0 1 0, 0 1 1, 1 1 1, 1 1 0, 0 1 0)),\
            ((0 0 1, 1 0 1, 1 1 1, 0 1 1, 0 0 1)))")
            .unwrap()
            .to_coordinates()
            .unwrap();
        let res = surface.as_x3d(0, false, None);
        assert_eq!(
            res.unwrap().geometry,
            "<IndexedFaceSet coordIndex='0 1 2 3 -1 4 5 6 7 -1 8 9 10 11 -1 12 13 14 15 -1 16 17 18 19 -1 20 21 22 23'>\
            <Coordinate point='0 0 0 0 0 1 0 1 1 0 1 0 0 0 0 0 1 0 1 1 0 1 0 0 0 0 0 1 0 0 1 0 1 0 0 1 1 1 0 1 1 1 1 0 1 \
            1 0 0 0 1 0 0 1 1 1 1 1 1 1 0 0 0 1 1 0 1 1 1 1 0 1 1' /></IndexedFaceSet>".to_string());
    }
}
