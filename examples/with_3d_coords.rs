extern crate sfcgal;
extern crate failure;

use sfcgal::{CoordSeq, ToCoordinates, ToSFCGAL};

fn fun() -> Result<(), failure::Error> {
    let coords_linestring = vec![(-0.5, -0.5, 2.5), (0., 0., 4.0)];
    let coords_polygon = vec![
        vec![(-1., -1., 3.0), (1., -1., 3.0), (1., 1., 3.0), (-1., 1., 3.0), (-1., -1., 3.0)], // Exterior ring
        vec![(0.1, 0.1, 3.0), (0.1, 0.9, 3.0), (0.9, 0.9, 3.0), (0.9, 0.1, 3.0), (0.1, 0.1, 3.0)], // 1 interior ring
    ];

    let line_3d = CoordSeq::Linestring(coords_linestring).to_sfcgal()?;
    let polygon_3d = CoordSeq::Polygon(coords_polygon).to_sfcgal()?;
    let intersects = line_3d.intersects_3d(&polygon_3d)?;
    assert!(intersects);
    let intersection = line_3d.intersection_3d(&polygon_3d)?;
    let coords_intersection: CoordSeq<(f64, f64, f64)> = intersection.to_coordinates()?;
    println!("{:?} and {:?} intersects at {:?}", line_3d, polygon_3d, coords_intersection);


    let coords_linestring: Vec<(f64, f64, f64)> = vec![];
    let coords_polygon = vec![
        vec![(-1., -1., 3.0), (1., -1., 3.0), (1., 1., 3.0), (-1., 1., 3.0), (-1., -1., 3.0)], // Exterior ring
        vec![(0.1, 0.1, 3.0), (0.1, 0.9, 3.0), (0.9, 0.9, 3.0), (0.9, 0.1, 3.0), (0.1, 0.1, 3.0)], // 1 interior ring
    ];

    let line_3d = CoordSeq::Linestring(coords_linestring).to_sfcgal()?;
    let polygon_3d = CoordSeq::Polygon(coords_polygon).to_sfcgal()?;
    let intersects = line_3d.intersects_3d(&polygon_3d)?;
    assert!(!intersects);
    let intersection = line_3d.intersection_3d(&polygon_3d)?;
    let _coords_intersection: CoordSeq<(f64, f64, f64)> = intersection.to_coordinates()?;
    println!("{:?} and {:?} intersects at {:?}", line_3d, polygon_3d, intersection);

    Ok(())
}

fn main() {
    fun().unwrap();
}
