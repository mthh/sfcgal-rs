use crate::{CoordSeq, Result};
use crate::coords::Point3d;

/// Formatted node element of a geometry in X3D xml
#[derive(Debug, PartialEq, Eq)]
pub struct X3dString(String);


/// Representation as [`X3dString`] (implemented on some ['CoordSeq`] variants
/// allowing to use it on some [`SFCGeometry`])
///
/// ['CoordSeq`]: struct.CoordSeq.html
/// [`SFCGeometry`]: struct.SFCGeometry.html
pub trait AsX3d {
    fn as_x3d(self, precision: usize, flip: bool, id: Option<String>) -> Result<X3dString>;
}

macro_rules! point_string {
    ($geom: expr, $precision: expr, $flip: expr) => ({
        if $flip == true {
            format!("{:.*} {:.*} {:.*}", $precision, $geom.1, $precision, $geom.0, $precision, $geom.2)
        } else {
            format!("{:.*} {:.*} {:.*}", $precision, $geom.0, $precision, $geom.1, $precision, $geom.2)
        }
    });
}

fn array_point_string(points: &[Point3d], precision: usize, flip: bool, skip_last: bool) -> String {
    points.iter()
        .enumerate()
        .filter_map(|(i, pt)| {
            if skip_last && i == points.len() -1 {
                None
            } else {
                Some(point_string!(pt, precision, flip))
            }
        })
        .collect::<Vec<String>>()
        .join(" ")
}

fn poly_string(rings: &[Vec<Point3d>], precision: usize, flip: bool) -> String {
    rings.iter()
        .map(|ring|{
            array_point_string(ring, precision, flip, true)
        })
        .collect::<Vec<String>>()
        .join(" ")
}


impl AsX3d for CoordSeq<Point3d> {
    fn as_x3d(self, precision: usize, flip: bool, id: Option<String>) -> Result<X3dString> {
        let id_obj = match id {
            None => String::default(),
            Some(_id) => format!("id='{}'", _id)
        };
        match self {
            // CoordSeq::Point(ref pt) => {
            //     Ok(X3dString(point_string!(pt, precision, flip, false)))
            // },
            // CoordSeq::Triangle(ref pts) => {
            //     Ok(X3dString(array_point_string(pts, precision, flip, false)))
            // },
            CoordSeq::Linestring(ref pts) => {
                let n_pts = pts.len();
                let s = format!(
                    "<LineSet {} vertexCount='{}'><Coordinate point='{}' /></LineSet>",
                    id_obj,
                    n_pts,
                    array_point_string(&pts, precision, flip, false),
                );
                Ok(X3dString(s))

            },
            CoordSeq::Triangulatedsurface(ref triangles) => {
                let n_triangles = triangles.len();
                let mut s = String::from(format!("<IndexedTriangleSet {} index='", id_obj));
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
                    s.push_str(&array_point_string(&triangles[i], precision, flip, false));
                    if i < n_triangles - 1 {
                        s.push_str(" ");
                    }
                }
                s.push_str("'/></IndexedTriangleSet>");
                Ok(X3dString(s))
            },
            CoordSeq::Polyhedralsurface(ref polygons) => {
                let n_geoms = polygons.len();
                let mut s = String::from(format!("<IndexedFaceSet {} coordIndex='", id_obj));
                let mut j = 0;
                polygons
                    .iter()
                    .enumerate()
                    .map(|(i, rings_patch)|{
                        let np: usize = rings_patch.iter().map(|g| g.len()).sum::<usize>() - 1usize;
                        for k in 0..np {
                            if (k > 0) { s.push_str(" "); }
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
                        poly_string(rings_patch.as_slice(), precision, flip)
                    })
                    .collect::<Vec<String>>()
                    .join(" "))
                );
                s.push_str("' /></IndexedFaceSet>");
                Ok(X3dString(s))
            },
            _ => Err(format_err!("Not implemented for type {:?}", self))
        }
    }
}



#[cfg(test)]
mod tests {
    use crate::{SFCGeometry, ToCoordinates, CoordSeq, ToSFCGAL, GeomType};
    use super::*;

    #[test]
    fn tin_representation() {
        let triangle = SFCGeometry::new("TIN (((0 0 0,0 0 1,0 1 0,0 0 0)), ((0 0 0,0 1 0,1 1 0,0 0 0)))")
            .unwrap()
            .to_coordinates()
            .unwrap();
        let res = triangle.as_x3d(0, false, Some("tin1".to_string()));
        assert_eq!(res.unwrap(), X3dString("<IndexedTriangleSet id='tin1' index='0 1 2 3 4 5'><Coordinate point='0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1 1 0'/></IndexedTriangleSet>".to_string()));
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
        let res = surface.as_x3d(0, false, Some("foo".to_string()));
        assert_eq!(res.unwrap(), X3dString("<IndexedFaceSet id='foo' coordIndex='0 1 2 3 -1 4 5 6 7 -1 8 9 10 11 -1 12 13 14 15 -1 16 17 18 19 -1 20 21 22 23'>\
<Coordinate point='0 0 0 0 0 1 0 1 1 0 1 0 0 0 0 0 1 0 1 1 0 1 0 0 0 0 0 1 0 0 1 0 1 0 0 1 1 1 0 1 1 1 1 0 1 1 0 0 0 1 0 0 1 1 1 1 1 1 1 0 0 0 1 1 0 1 1 1 1 0 1 1' />\
</IndexedFaceSet>".to_string()));
    }
}
