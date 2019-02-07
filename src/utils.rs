use std::ffi::CStr;
use libc::c_char;
use sfcgal_sys::sfcgal_geometry_t;
use approx::AbsDiff;

use crate::Result;
use crate::errors::get_last_error;

pub(crate) fn _string(raw_ptr: *const c_char) -> String {
    let c_str = unsafe { CStr::from_ptr(raw_ptr) };
    std::str::from_utf8(c_str.to_bytes()).unwrap().to_string()
}

pub(crate) fn check_null_geom(g: *const sfcgal_geometry_t) -> Result<()> {
    if g.is_null() {
        return Err(format_err!(
             "Error - Encoutered a null Geometry : {}",
             get_last_error()
        ));
    }
    Ok(())
}

pub(crate) fn check_predicate(val: i32) -> Result<bool> {
    match val {
        1 => Ok(true),
        0 => Ok(false),
        _ => Err(format_err!("SFCGAL error: {}", get_last_error()))
    }
}

pub(crate) fn check_computed_value(val: f64) -> Result<f64> {
    match AbsDiff::default().eq(&val, &-1.0) {
        true => Err(format_err!("SFCGAL error: {}", get_last_error())),
        false => Ok(val)
    }
}
