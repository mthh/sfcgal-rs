# Changelog
All notable changes to this project will be documented in this file.

## Unreleased
### Added
- Update to SFCGAL 1.4.1 and expose the new following functions from the C API:
  * `sfcgal_geometry_alpha_shapes` and `sfcgal_geometry_optimal_alpha_shapes`.

### Changed
- Update various dependencies (notably `geo_types` and `geojson`).
- Switch from travis to github actions.
- Apply some clippy suggestions.


## [0.5.0] - 2020-10-18
### Added
- Update to SFCGAL 1.3.8 and expose the following functions from the C API:
  * `sfcgal_geometry_covers`, `sfcgal_geometry_covers_3d`, `sfcgal_geometry_is_measured` and `sfcgal_geometry_line_sub_string`


## [0.4.0] - 2020-09-02
### Added
- Enable conversion from/to geo_types Triangle and conversion from geo_types Rect

### Changed
- Update `geo` (to 0.6) and `sfcgal-sys` (to 0.4).


## [0.3.0] - 2020-07-08
### Changed
- Use ffi types : size_t, c_char, etc [#1](https://github.com/mthh/sfcgal-rs/pull/1)
- Use less restrictive x.y version number in Cargo.toml [#1](https://github.com/mthh/sfcgal-rs/pull/1)
- Replace deprecated uninitialized() with MaybeUninit [#1](https://github.com/mthh/sfcgal-rs/pull/1)
- Tidy "use ..." statements and other code simplifications [#1](https://github.com/mthh/sfcgal-rs/pull/1)
- Replace deprecated `failure` crate with `anyhow` [#1](https://github.com/mthh/sfcgal-rs/pull/1)
- Update `geo` (to 0.5), `geojson` (0.19) and `sfcgal-sys` (to 0.3).
