extern crate azul;
extern crate svg as osvg;

extern crate geojson;
extern crate geo_types;
extern crate sfcgal;

use geojson::{GeoJson, Geometry, Feature, FeatureCollection, Value, conversion::TryInto as TryIntoGeoType};
use sfcgal::{ToSfcgal, TryInto};
use std::fs::File;
use std::io::Read;

use azul::{
    prelude::*,
    widgets::{button::Button, svg::*},
    dialogs::*,
};

use osvg::Document;
use osvg::Node;
use osvg::node::element::{Group, Path};
use osvg::node::element::path::Data;

#[derive(Debug)]
pub struct MyAppData {
    pub svg: Option<(SvgCache<MyAppData>, Vec<SvgLayerId>)>,
}

impl Layout for MyAppData {
    fn layout(&self, info: WindowInfo)
    -> Dom<MyAppData>
    {
        if let Some((svg_cache, svg_layers)) = &self.svg {
            Svg::with_layers(svg_layers).dom(&info.window, &svg_cache)
        } else {
            Dom::new(NodeType::Div)
                .with_class("__azul-native-button")
                .with_callback(On::MouseUp, Callback(my_button_click_handler))
        }
    }
}

fn my_button_click_handler(app_state: &mut AppState<MyAppData>, _event: WindowEvent) -> UpdateScreen {
    const TEST_SVG: &str = include_str!("export.svg");
    let mut svg_cache = SvgCache::empty();
    let svg_layers = svg_cache.add_svg(TEST_SVG).unwrap();
    app_state.data.modify(|data| data.svg = Some((svg_cache, svg_layers)));
    UpdateScreen::Redraw
}


fn main() {
    let (bg, res) = make_skeleton();
    create_svg_from_geojson(bg, res);
    let app = App::new(MyAppData { svg: None }, AppConfig::default());
    app.run(Window::new(WindowCreateOptions::default(), css::native()).unwrap()).unwrap();
}

fn make_skeleton() -> (GeoJson, GeoJson) {
    let path = "examples/abc.geojson";
    let mut file = File::open(path).unwrap_or_else(|err| {
        println!("Unable to open layer at path: \"{}\"\nError: {}", path, err);
        std::process::exit(1)
    });
    let mut raw_json = String::new();
    file.read_to_string(&mut raw_json).unwrap();
    let decoded_geojson = raw_json.parse::<GeoJson>().unwrap();
    let features = match decoded_geojson {
        GeoJson::FeatureCollection(ref collection) => &collection.features,
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
    });
    (decoded_geojson, result)
}

pub fn create_svg_from_geojson(background_geojson: GeoJson, result: GeoJson) {
    let mapextent = MapExtent {
      left: -80.,
      right: -72.,
      bottom: 4.,
      top: 9.,
    };
    let width: u32 = 600;
    let height: u32 = 600;
    let converter = Converter::new(width, height, &mapextent);
    // Create a new svg document:
    let mut document = Document::new()
        .set("x", "0")
        .set("y", "0")
        .set("width", format!("{}", converter.viewport_width))
        .set("height", format!("{}", converter.viewport_height));

    let features = match background_geojson {
        GeoJson::FeatureCollection(collection) => collection.features,
        _ => {
            println!("Expected a GeoJSON feature collection!");
            std::process::exit(1)
        }
    };

    let mut group = Group::new();
    for feature in features {
        let geom = feature.geometry.unwrap();
        match geom.value {
            Value::Polygon(positions) => {
                group.append(Path::new()
                    .set("fill", "gray")
                    .set("fill-opacity", "0.8")
                    .set("stroke", "blue")
                    .set("stroke-width", "0.8")
                    .set("stroke-opacity", "0.9")
                    .set("d", converter.draw_path_ring(&positions, None)));
            }
            Value::MultiPolygon(polys) => {
                let mut data = Data::new();
                for positions in &polys {
                    data = converter.draw_path_ring(positions, Some(data));
                }
                data = data.close();
                group.append(Path::new()
                    .set("fill", "gray")
                    .set("fill-opacity", "0.8")
                    .set("stroke", "blue")
                    .set("stroke-width", "0.8")
                    .set("stroke-opacity", "0.9")
                    .set("d", data));
            }
            _ => panic!("Error: Expected a polygon / multipolygon!!"),
        }
    }
    document = document.add(group.set("id", "background"));

    let features = match result {
        GeoJson::FeatureCollection(collection) => collection.features,
        _ => {
            println!("Expected a GeoJSON feature collection!");
            std::process::exit(1)
        }
    };

    let mut group = Group::new();
    for feature in features {
        let geom = feature.geometry.unwrap();
        match geom.value {
            Value::MultiLineString(ref lines) => {
                let mut data = Data::new();
                for positions in lines {
                    data = converter.draw_path_ring(&[positions.to_vec()], Some(data));
                }
                group.append(Path::new()
                                 .set("fill", "none")
                                 .set("stroke", "orange")
                                 .set("stroke-width", "1.1")
                                 .set("stroke-opacity", "0.9")
                                 .set("d", data));
            }
            _ => panic!("Error: Expected a linestring / multilinestring!!"),
        }
    }
    document = document.add(group.set("id", "background"));
    svg::save("examples/export.svg", &document).unwrap();
}


#[derive(Debug, Default, Clone)]
pub struct MapExtent {
    pub left: f64,
    pub right: f64,
    pub bottom: f64,
    pub top: f64,
}

struct Converter<'a> {
    viewport_width: u32,
    viewport_height: u32,
    map_extent: &'a MapExtent,
    resolution: f64,
}

impl<'a> Converter<'a> {
    pub fn new(viewport_width: u32, viewport_height: u32, map_extent: &'a MapExtent) -> Self {
        let xres = (map_extent.right - map_extent.left) / viewport_width as f64;
        let yres = (map_extent.top - map_extent.bottom) / viewport_height as f64;
        let res = xres.max(yres);
        Converter {
            viewport_width: viewport_width,
            viewport_height: viewport_height,
            map_extent: map_extent,
            resolution: res,
        }
    }

    pub fn draw_path_ring(&self, positions: &[Vec<Vec<f64>>], d: Option<Data>) -> Data {
        let (mut data, close) = if d.is_some() {
            (d.unwrap(), false)
        } else {
            (Data::new(), true)
        };
        for ring in positions {
            let mut iter = ring.iter();
            let first = iter.next().unwrap();
            data = data.move_to(((first[0] - self.map_extent.left) / self.resolution,
                                 (self.map_extent.top - first[1]) / self.resolution));
            for point in iter {
                data = data.line_to(((point[0] - self.map_extent.left) / self.resolution,
                                     (self.map_extent.top - point[1]) / self.resolution));
            }
        }
        if close { data.close() } else { data }
    }
}
