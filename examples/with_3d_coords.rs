extern crate sfcgal;
extern crate failure;

use sfcgal::{CoordSeq, ToCoordinates, ToSFCGAL};

fn fun() -> Result<(), failure::Error> {
    let coords_linestring = vec![(-0.5, -0.5, 2.5), (0., 0., 4.0)];
    let coords_polygon = vec![
        vec![(-1., -1., 3.0), (1., -1., 3.0), (1., 1., 3.0), (-1., 1., 3.0), (-1., -1., 3.0)], // Exterior ring
        vec![(0.1, 0.1, 3.0), (0.1, 0.9, 3.0), (0.9, 0.9, 3.0), (0.9, 0.1, 3.0), (0.1, 0.1, 3.0)], // 1 interior ring
    ];

    let line = CoordSeq::Linestring(coords_linestring).to_sfcgal()?;
    let poly = CoordSeq::Polygon(coords_polygon).to_sfcgal()?;
    let intersects = line.intersects_3d(&poly)?;
    assert!(intersects);
    let intersection = line.intersection_3d(&poly)?;
    let coordinates_intersection: CoordSeq<(f64, f64, f64)> = intersection.to_coordinates()?;
    println!("{:?} and {:?} interescts at coordinates {:?}", line, poly, coordinates_intersection);
    Ok(())
}

fn main() {
    fun().unwrap();
}
