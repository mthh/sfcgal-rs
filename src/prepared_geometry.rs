// use sfcgal_sys::{
//     initialize, sfcgal_prepared_geometry_t,
//     sfcgal_prepared_geometry_create, sfcgal_prepared_geometry_delete,
// };

// #[repr(C)]
// pub struct SFCPreparedGeometry(NonNull<sfcgal_prepared_geometry_t>);
//
// impl Drop for SFCPreparedGeometry {
//     fn drop(&mut self) {
//         unsafe { sfcgal_prepared_geometry_delete(self.0.as_mut()) }
//     }
// }
