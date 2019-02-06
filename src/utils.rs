use std::ffi::CStr;
use libc::c_char;

pub fn _string(raw_ptr: *const c_char) -> String {
    let c_str = unsafe { CStr::from_ptr(raw_ptr) };
    std::str::from_utf8(c_str.to_bytes()).unwrap().to_string()
}
