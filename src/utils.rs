use std::{ffi::CStr, os::raw::c_char};

use approx::AbsDiff;
use sfcgal_sys::{sfcgal_geometry_t, sfcgal_prepared_geometry_t};

use crate::errors::get_last_error;
use crate::Result;

pub(crate) fn _string(raw_ptr: *const c_char) -> String {
    let c_str = unsafe { CStr::from_ptr(raw_ptr) };

    std::str::from_utf8(c_str.to_bytes()).unwrap().to_string()
}

pub(crate) fn check_null_geom(g: *const sfcgal_geometry_t) -> Result<()> {
    if g.is_null() {
        return Err(format_err!(
            "Error - Encountered a null Geometry : {}",
            get_last_error()
        ));
    }

    Ok(())
}

pub(crate) fn check_null_prepared_geom(g: *mut sfcgal_prepared_geometry_t) -> Result<()> {
    if g.is_null() {
        return Err(format_err!(
            "Error - Encountered a null Prepared Geometry : {}",
            get_last_error()
        ));
    }

    Ok(())
}

pub(crate) fn check_predicate(val: i32) -> Result<bool> {
    match val {
        1 => Ok(true),
        0 => Ok(false),
        _ => Err(format_err!("SFCGAL error: {}", get_last_error())),
    }
}

pub(crate) fn check_computed_value(val: f64) -> Result<f64> {
    if AbsDiff::default().eq(&val, &-1.0) {
        Err(format_err!("SFCGAL error: {}", get_last_error()))
    } else {
        Ok(val)
    }
}

pub(crate) fn check_nan_value(val: f64) -> Result<f64> {
    if val.is_nan() {
        Err(format_err!("SFCGAL error: {}", get_last_error()))
    } else {
        Ok(val)
    }
}

// Use the size returned by the C API to build the string
// (as it seems to not always end with a null byte)
// from the pointer to uninitialized memory with give
// to it earlier.
pub(crate) fn _c_string_with_size(raw_ptr: *const c_char, size: usize) -> String {
    let slice: &[u8] = unsafe { std::slice::from_raw_parts(raw_ptr as *const u8, size) };

    let res = std::str::from_utf8(slice).unwrap().to_string();

    unsafe { libc::free(raw_ptr as *mut libc::c_void) };

    res
}

// Get the raw u8 bytes
// The users can proceed to use it however they want.
pub(crate) fn get_raw_bytes(raw_ptr: *const c_char, size: usize) -> Vec<u8> {
    let slice: &[u8] = unsafe { std::slice::from_raw_parts(raw_ptr as *const u8, size) };

    slice.to_vec()
}
