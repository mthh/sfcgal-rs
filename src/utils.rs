use std::ffi::CStr;
use libc::c_char;
use sfcgal_sys::sfcgal_geometry_t;

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
