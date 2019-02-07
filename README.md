# sfcgal-rs
__*(WIP)*__

[![Build Status Travis](https://travis-ci.org/mthh/sfcgal-rs.svg?branch=master)](https://travis-ci.org/mthh/sfcgal-rs)

Rust bindings to [`SFCGAL`](http://oslandia.github.io/SFCGAL/) C API.  
Based on the [sfcgal-sys](https://github.com/mthh/sfcgal-rs) crate exposing low-level bindings.

## Features / TODO

- [x] `sfcgal_geometry_*`  
- [x] Conversion from / to [geo-types](https://github.com/georust/geo) with working examples
- [ ] `sfcgal_prepared_geometry_*`  
- [ ] Nice documentation
- [ ] Warning messages

## Usage
```rust
extern crate geo_types;
extern crate sfcgal;

use geo_types::{LineString, Polygon};
use sfcgal::ToSfcgal;


// create a geo_types Polygon:
let polygon = Polygon::new(
    LineString::from(vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.,), (0., 0.)]),
    vec![LineString::from(
        vec![(0.1, 0.1), (0.1, 0.9,), (0.9, 0.9), (0.9, 0.1), (0.1, 0.1)])]);

// create a geo_types LineString:
let line = LineString::from(vec![(-0.5, -0.5), (1.3, 1.), (1.1, 1.9)]);

// convert them to sfcgal geometries:
let polyg_sfc = polygon.to_sfcgal().unwrap();
let line_sfc = line.to_sfcgal().unwrap();

// Use SFCGAL operations:
assert!(polyg_sfc.intersects(&line_sfc).unwrap(), true);
```

### Examples

See `examples/skeleton_geojson.rs` for an example of working with some other crates from Rust geo ecosystem.  
See `examples/render_skeleton.rs` for an example of computing a straight skeleton from some GeoJSON file and rendering the result with [web_view](https://github.com/Boscop/web-view/) after converting it in SVG.  


## License

Licensed under either of
 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you shall be dual licensed as above, without any
additional terms or conditions.
