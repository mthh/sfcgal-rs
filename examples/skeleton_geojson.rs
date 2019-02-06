extern crate geojson;
extern crate geo_types;
extern crate sfcgal;

use geojson::{GeoJson, Geometry, Feature, FeatureCollection, Value, conversion::TryInto as TryIntoGeoType};
use sfcgal::{ToSfcgal, TryInto};
use std::fs::File;
use std::io::Read;

fn main() {
    let path = "examples/abc.geojson";
    let mut file = File::open(path).unwrap_or_else(|err| {
        println!("Unable to open layer at path: \"{}\"\nError: {}", path, err);
        std::process::exit(1)
    });
    let mut raw_json = String::new();
    file.read_to_string(&mut raw_json).unwrap();
    let decoded_geojson = raw_json.parse::<GeoJson>().unwrap();
    let features = match decoded_geojson {
        GeoJson::FeatureCollection(collection) => collection.features,
        _ => panic!("Error: expected a Feature collection!"),
    };
    let mut sfcgal_geoms = vec![];
    features
        .iter()
        .map(|ref feature| {
            if let Some(geom) = &feature.geometry {
                if let Value::Polygon(ref pos) = geom.value {
                    // let geojson_polygon = Value::Polygon(positions);
                    let geo_polygon: geo_types::Polygon<f64> = Value::Polygon(pos.to_vec()).try_into().unwrap();
                    sfcgal_geoms.push(geo_polygon.to_sfcgal().unwrap());
                }
            }
        }).for_each(drop);

    let results = sfcgal_geoms
        .iter()
        .map(|sfc_geom| {
            sfc_geom.straight_skeleton().unwrap()
        }).collect::<Vec<_>>();

    let features_geojson = results
        .into_iter()
        .map(|res_sfc_geom|{
            let geogeom: geo_types::MultiLineString<f64> = res_sfc_geom
                .try_into()
                .unwrap()
                .as_multilinestring()
                .unwrap();
            Feature {
                bbox: None,
                geometry: Some(Geometry {
                    value: Value::from(&geogeom),
                    bbox: None,
                    foreign_members: None,
                }),
                id: None,
                properties: None,
                foreign_members: None,
            }
        }).collect::<Vec<_>>();

    let result = GeoJson::FeatureCollection(FeatureCollection {
        foreign_members: None,
        bbox: None,
        features: features_geojson,
    }).to_string();
    print!("{}", result);
}
