# sfcgal-rs

[![Build and run tests](https://github.com/mthh/sfcgal-rs/actions/workflows/test.yml/badge.svg)](https://github.com/mthh/sfcgal-rs/actions/workflows/test.yml)
[![Crates.io](https://img.shields.io/crates/v/sfcgal.svg)](https://crates.io/crates/sfcgal)
[![Documentation](https://img.shields.io/badge/documentation-0.8.2-green)](https://mthh.github.io/sfcgal-rs/sfcgal/)


Rust bindings providing a high-level API to [`SFCGAL`](https://sfcgal.gitlab.io/SFCGAL/) library and conversion to / from other geometry crates from Rust ecosystem.  
Based on the [sfcgal-sys](https://github.com/mthh/sfcgal-sys) crate exposing low-level bindings.

Some key features of the underlying library:
- Supports ISO 19107 and [OGC Simple Features Access 1.2](http://www.opengeospatial.org/standards/sfa) for 3D operations.
- Reads and writes WKT with exact rational number representation of coordinates for 2D and 3D geometries.
- Intersection, difference and union.
- Straight skeleton, tesselation, Minkovski sum, alpha shapes and convex hull.

> [!IMPORTANT]
> Note that the required version of SFCGAL is currently 2.x.x (latest version on 2025-05-27 is 2.1.0).  
> If you want to use SFCCAL 1.5.x, you can use the 0.7.x version of this crate.  
> If you want to use SFCCAL 1.4.x, you can use the 0.6.x version of this crate.

## Usage

__Example with 3-member tuples for 3d coordinates and WKT__:
```rust
extern crate sfcgal;
use sfcgal::{SFCGeometry, CoordSeq, ToCoordinates, ToSFCGAL};

// create a linestring from WKT:
let line_3d = SFCGeometry::new("LINESTRING(-0.5 -0.5 2.5, 0.0 0.0 4.0)")?;

// create a polygon as Vec of 3-member tuples...
let coords_polygon = vec![
    vec![ // Exterior ring
        (-1., -1., 3.0),
        (1., -1., 3.0),
        (1., 1., 3.0),
        (-1., 1., 3.0),
        (-1., -1., 3.0),
    ],
    vec![ // 1 interior ring
        (0.1, 0.1, 3.0),
        (0.1, 0.9, 3.0),
        (0.9, 0.9, 3.0),
        (0.9, 0.1, 3.0),
        (0.1, 0.1, 3.0),
    ],
];
// ...by using the CoordSeq enum variants to match the wanted SFCGAL geometry type
// (returns a SFCGeometry)
let polygon_3d = CoordSeq::Polygon(coords_polygon).to_sfcgal()?;

// ...
let intersection = line_3d.intersection_3d(&polygon_3d)?;

// Retrieve coordinates of the resulting geometry as 3-member tuples:
let coords_intersection: CoordSeq<(f64, f64, f64)> = intersection.to_coordinates()?;

println!("{:?} and {:?} intersects at {:?}", line_3d, polygon_3d, coords_intersection);
```

__Example with [geo-types](https://github.com/georust/geo)__:
```rust
extern crate geo_types;
extern crate sfcgal;

use geo_types::{LineString, Polygon};
use sfcgal::ToSFCGAL;


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

See `examples/skeleton_geojson.rs` for an example of working with some other crates from Rust-geo ecosystem.  


### Motivation

Needed a SFCGAL feature for a side-project in Rust and I thought it would be a good opportunity to try using [bindgen](https://github.com/rust-lang/rust-bindgen) on SFCGAL C API.  
In the end a large part of the API is wrapped, so maybe it could be reused or improved by the community now it's published on [crates.io](https://crates.io/crates/sfcgal).

## Other SFCGAL bindings

Those Rust bindings are community-driven. If you are looking for official bindings, you can check the [PySFCGAL](https://gitlab.com/sfcgal/pysfcgal) project
which is developed by the SFCGAL team.

## License

Licensed under either of
 * Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you shall be dual licensed as above, without any
additional terms or conditions.
