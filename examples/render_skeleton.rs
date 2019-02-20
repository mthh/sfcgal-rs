#![windows_subsystem = "windows"]

extern crate web_view;
extern crate geojson;
extern crate geo_types;
extern crate sfcgal;
extern crate svg;

use web_view::*;
use geojson::{GeoJson, Geometry, Feature, FeatureCollection, Value, conversion::TryInto as TryIntoGeoType};
use sfcgal::{ToSFCGAL, TryInto};
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::thread;

use svg::Document;
use svg::Node;
use svg::node::element::{Group, Path};
use svg::node::element::path::Data;


fn read_example_file() -> String {
    let path = "examples/abc.geojson";
    let mut file = File::open(path).unwrap_or_else(|err| {
        println!("Unable to open layer at path: \"{}\"\nError: {}", path, err);
        std::process::exit(1)
    });
    let mut raw_json = String::new();
    file.read_to_string(&mut raw_json).unwrap();
    raw_json
}

fn make_skeleton(raw_json: String) -> (GeoJson, GeoJson) {
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
                    let geo_polygon: geo_types::Polygon<f64> = Value::Polygon(pos.to_vec()).try_into().unwrap();
                    sfcgal_geoms.push(geo_polygon.to_sfcgal().unwrap());
                }
            }
        }).for_each(drop);

    let features_geojson = sfcgal_geoms
        .into_iter()
        .map(|sfc_geom| {
            let res_sfc_geom = sfc_geom.straight_skeleton().unwrap();
            let geogeom: geo_types::MultiLineString<f64> = res_sfc_geom
                .try_into()
                .unwrap()
                .into_multi_line_string()
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

    let res = GeoJson::FeatureCollection(FeatureCollection {
        foreign_members: None,
        bbox: None,
        features: features_geojson,
    });
    (decoded_geojson, res)
}

pub fn create_svg_from_geojson(layers: &[(GeoJson, String)], _width: u32, _height: u32) -> String {
    let mapextent = MapExtent {
      left: -78.5,
      right: -73.,
      bottom: 5.,
      top: 9.,
    };
    let width: u32 = _width - 20;
    let height: u32 = _height - 20;
    let converter = Converter::new(width, height, mapextent);
    // Create a new svg document:
    let mut document = Document::new()
        .set("x", "0")
        .set("y", "0")
        .set("width", format!("{}", converter.viewport_width))
        .set("height", format!("{}", converter.viewport_height));

    let mut result_string = String::new();

    for (geojson_lyr, id_layer) in layers {
        let features = match geojson_lyr {
            GeoJson::FeatureCollection(collection) => &collection.features,
            _ => {
                println!("Expected a GeoJSON feature collection!");
                std::process::exit(1)
            }
        };
        let mut group = Group::new();
        for feature in features {
            if let Some(ref geom) = feature.geometry {
                match &geom.value {
                    Value::Polygon(ref positions) => {
                        group.append(Path::new()
                            .set("fill", "gray")
                            .set("fill-opacity", "0.8")
                            .set("stroke", "blue")
                            .set("stroke-width", "0.8")
                            .set("stroke-opacity", "0.9")
                            .set("d", converter.draw_path_ring(&positions, None)));
                    }
                    Value::MultiPolygon(ref polys) => {
                        let mut data = Data::new();
                        for positions in polys {
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
                    _ => panic!("Error: Expected a multilinestring / polygon / multipolygon!!"),
                }
            }
        }
        document = document.add(group.set("id", id_layer.as_str()));
    }
    svg::write(unsafe { BufWriter::new(result_string.as_mut_vec()) }, &document).unwrap();
    result_string
}


#[derive(Debug, Default, Clone)]
pub struct MapExtent {
    pub left: f64,
    pub right: f64,
    pub bottom: f64,
    pub top: f64,
}

struct Converter {
    viewport_width: u32,
    viewport_height: u32,
    map_extent: MapExtent,
    resolution: f64,
}

impl Converter {
    pub fn new(viewport_width: u32, viewport_height: u32, map_extent: MapExtent) -> Self {
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

fn main() {
    let raw_json = read_example_file();
    let decoded_geojson = raw_json.parse::<GeoJson>().unwrap();
    let webview = web_view::builder()
        .title("Straight skeleton example")
        .content(Content::Html(HTML))
        .size(800, 700)
        .resizable(true)
        .debug(true)
        .user_data(0)
        .invoke_handler(|webview, arg| {
            match arg {
                "draw" => {
                    webview.eval("document.getElementById('waiting').style.display = '';").unwrap();
                    let handle = webview.handle();
                    let raw_json_clone = raw_json.clone();
                    thread::spawn(move || {
                        let (decoded_geojson, res) = make_skeleton(raw_json_clone);
                        let svg_content = create_svg_from_geojson(
                            &vec![
                                (decoded_geojson, String::from("bg")),
                                (res, String::from("result")),
                            ], 800, 600);
                        handle.dispatch(move |webview| {
                            render(webview, svg_content)
                        }).unwrap();
                    });
                }
                "exit" => {
                    webview.terminate();
                }
                _ => unimplemented!(),
            };
            Ok(())
        })
        .build()
        .unwrap();

    webview.handle().dispatch(move |webview| {
        render(webview, create_svg_from_geojson(
            &vec![(decoded_geojson, String::from("bg"))], 800, 600))
    }).unwrap();

    webview.run().unwrap();
}

fn render(webview: &mut WebView<i32>, svg_content: String) -> WVResult {
    webview.eval(&format!("update_svg(`{}`)", svg_content))
}

const HTML: &str = r#"
<!doctype html>
<html>
    <body>
        <div style="margin-left: 20px; margin-right: 20px;">
        	<button style="float: left; margin: 20px;" onclick="external.invoke('draw')">Draw</button>
        	<button style="float: right; margin: 20px;" onclick="external.invoke('exit')">Exit</button>
        </div>
        <div style="text-align: center; clear: both;">
            <p style="margin: auto; display: none;" id="waiting">Working...</p>
        </div>
    	<div id="result"></div>
    	<script type="text/javascript">
    		function update_svg(n) {
    			document.getElementById('result').innerHTML = n;
                document.getElementById('waiting').style.display = 'none';
    		}
    	</script>
    </body>
</html>
"#;
