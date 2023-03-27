use crate::utils::_string;
use anyhow::Error;
use sfcgal_sys::w_sfcgal_get_last_error;

/// Type alias for [`Result`] with the error type set to [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

pub fn get_last_error() -> String {
    let message = unsafe { w_sfcgal_get_last_error() };
    _string(message)
}
