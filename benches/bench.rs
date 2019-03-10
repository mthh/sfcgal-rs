#![feature(test)]
extern crate geo_types;
extern crate geojson;
extern crate sfcgal;
extern crate test;

use geojson::conversion::TryInto as TryIntoGeoType;
use sfcgal::{FromGeoJSON, SFCGeometry, ToCoordinates, ToGeoJSON, ToSFCGAL, TryInto};
use std::io::Read;
use test::{black_box, Bencher};

type Point2d = (f64, f64);
type Point3d = (f64, f64, f64);

fn read_example_file() -> Vec<geojson::Geometry> {
    let path = "examples/abc.geojson";
    let mut file = std::fs::File::open(path).unwrap_or_else(|err| {
        println!("Unable to open layer at path: \"{}\"\nError: {}", path, err);
        std::process::exit(1)
    });
    let mut raw_json = String::new();
    file.read_to_string(&mut raw_json).unwrap();
    let decoded_geojson = raw_json.parse::<geojson::GeoJson>().unwrap();
    let features = match decoded_geojson {
        geojson::GeoJson::FeatureCollection(ref collection) => &collection.features,
        _ => panic!("Error: expected a Feature collection!"),
    };
    let mut geoms = Vec::new();
    features
        .iter()
        .map(|feature| {
            if feature.geometry.is_some() {
                geoms.push(feature.geometry.clone().take().unwrap());
            }
        })
        .for_each(drop);
    geoms
}

fn make_sfgal_geom() -> SFCGeometry {
    let geoms = read_example_file();
    SFCGeometry::from_geojson::<Point2d>(&geoms[0].value).unwrap()
}

#[bench]
fn bench_geojson_to_geotypes_to_sfcgal_2d(b: &mut Bencher) {
    let mut geom = read_example_file();
    let g = geom.get_mut(0).unwrap();
    b.iter(|| {
        let geo_polygon: geo_types::Polygon<f64> = g.value.clone().try_into().unwrap();
        geo_polygon.to_sfcgal().unwrap();
    });
}

#[bench]
fn bench_geotypes_to_sfcgal_2d(b: &mut Bencher) {
    let mut geom = read_example_file();
    let g = geom.get_mut(0).unwrap();
    let geo_polygon: geo_types::Polygon<f64> = g.value.clone().try_into().unwrap();
    b.iter(|| {
        geo_polygon.to_sfcgal().unwrap();
    });
}

#[bench]
fn bench_geojson_to_sfcgal_3d(b: &mut Bencher) {
    let geom = read_example_file();
    b.iter(|| SFCGeometry::from_geojson::<Point3d>(&geom[0].value).unwrap());
}

#[bench]
fn bench_geojson_to_sfcgal_2d(b: &mut Bencher) {
    let geom = read_example_file();
    b.iter(|| SFCGeometry::from_geojson::<Point2d>(&geom[0].value).unwrap());
}

#[bench]
fn bench_sfcgal_2d_to_coordinates(b: &mut Bencher) {
    let geom = make_sfgal_geom();
    b.iter(|| geom.to_coordinates::<Point2d>().unwrap());
}

#[bench]
fn bench_sfcgal_2d_to_geojson(b: &mut Bencher) {
    let geom = make_sfgal_geom();
    b.iter(|| {
        let _coords = geom.to_geojson::<Point2d>().unwrap();
    });
}

#[bench]
fn bench_sfcgal_2d_to_geotype_geometry(b: &mut Bencher) {
    let geom = make_sfgal_geom();
    b.iter(|| {
        let _polyg: geo_types::Geometry<f64> = geom.clone().try_into().unwrap();
    });
}

#[bench]
fn bench_sfcgal_2d_to_geotype_concrete_type(b: &mut Bencher) {
    let geom = make_sfgal_geom();
    b.iter(|| {
        let _polyg: geo_types::Polygon<f64> =
            geom.clone().try_into().unwrap().into_polygon().unwrap();
    });
}
