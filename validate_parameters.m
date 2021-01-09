function [point_ids, calib_points, global_r_prctile] = validate_parameters(outdir, data, dim_mfld, se_targ, params)
    
    validate_outdir(outdir)
    validate_data(data, dim_mfld);
    validate_se(se_targ);
  
    point_ids    = sort(params.Results.point_ids);
    if (ismember('calib_points', params.UsingDefaults))
        calib_points = point_ids;
    else
        calib_points = sort(params.Results.calib_points);
    end
    global_r_prctile = params.Results.global_r_prctile;
    
end
