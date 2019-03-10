extern crate sfcgal;
use sfcgal::{SFCGeometry, CoordSeq, ToCoordinates, ToSFCGAL};

fn main() {
    // create a linestring from WKT:
    let line_3d = SFCGeometry::new("LINESTRING(-0.5 -0.5 2.5, 0.0 0.0 4.0)").unwrap();

    // create a polygon as Vec of 3-member tuples...
    let coords_polygon = vec![
        vec![(-1., -1., 3.0), (1., -1., 3.0), (1., 1., 3.0), (-1., 1., 3.0), (-1., -1., 3.0)], // Exterior ring
        vec![(0.1, 0.1, 3.0), (0.1, 0.9, 3.0), (0.9, 0.9, 3.0), (0.9, 0.1, 3.0), (0.1, 0.1, 3.0)], // 1 interior ring
    ];
    // ...by using the CoordSeq enum variants to match the wanted SFCGAL geometry type
    // (returns a SFCGeometry)
    let polygon_3d = CoordSeq::Polygon(coords_polygon).to_sfcgal().unwrap();

    let intersects = line_3d.intersects_3d(&polygon_3d).unwrap();
    assert!(intersects);
    let intersection = line_3d.intersection_3d(&polygon_3d).unwrap();
    let coords_intersection: CoordSeq<(f64, f64, f64)> = intersection.to_coordinates().unwrap();
    println!("{:?} and {:?} intersects at {:?}", line_3d, polygon_3d, coords_intersection);
}
