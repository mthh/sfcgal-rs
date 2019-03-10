extern crate geo_types;
extern crate geojson;
extern crate sfcgal;

use geojson::{Feature, FeatureCollection, GeoJson, Geometry, Value};
use sfcgal::{FromGeoJSON, SFCGeometry, ToGeoJSON};
use std::{fs::File, io::Read};

type Pt3d = (f64, f64, f64);

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
    let new_features = features
        .iter()
        .filter_map(|ref feature| {
            if let Some(geom) = &feature.geometry {
                if let Value::Polygon(..) = geom.value {
                    let polygon_sfc = SFCGeometry::from_geojson::<Pt3d>(&geom.value).unwrap();
                    let skeleton = polygon_sfc.straight_skeleton().unwrap();
                    let geojson_geom: Value = skeleton.to_geojson::<Pt3d>().unwrap();
                    return Some(Feature {
                        bbox: None,
                        geometry: Some(Geometry {
                            value: geojson_geom,
                            bbox: None,
                            foreign_members: None,
                        }),
                        id: None,
                        properties: feature.properties.clone(),
                        foreign_members: None,
                    });
                }
            }
            return None;
        })
        .collect::<Vec<_>>();

    let _result_geojson = GeoJson::FeatureCollection(FeatureCollection {
        foreign_members: None,
        bbox: None,
        features: new_features,
    })
    .to_string();
}
