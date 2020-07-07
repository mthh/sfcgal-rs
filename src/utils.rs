use approx::AbsDiff;
use libc::c_char;
use sfcgal_sys::sfcgal_geometry_t;
use std::ffi::CStr;

use crate::errors::get_last_error;
use crate::Result;

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

// Use the size returned by the C API to build the string
// (as it seems to not always end with a null byte)
// from the pointer to uninitialized memory with give
// to it earlier.
pub(crate) fn _c_string_with_size(raw_ptr: *mut c_char, size: usize) -> String {
    let slice: &[u8] =
        unsafe { &*(std::slice::from_raw_parts(raw_ptr, size) as *const [i8] as *const [u8]) };
    let res = std::str::from_utf8(slice).unwrap().to_string();
    unsafe { libc::free(raw_ptr as *mut libc::c_void) };
    res
}
