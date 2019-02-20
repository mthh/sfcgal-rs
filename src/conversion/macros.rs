macro_rules! make_sfcgal_multi_geom {
    ($c_geom: expr, $iter: expr) => ({
        let out_multi = unsafe { $c_geom };
        check_null_geom(out_multi)?;
        for single_sfcgal_geom in $iter {
            unsafe {
                sfcgal_geometry_collection_add_geometry(out_multi, single_sfcgal_geom as *mut sfcgal_geometry_t)
            };
        }
        unsafe { SFCGeometry::new_from_raw(out_multi, true) }
    });
}
