function [point_ids, calib_points, ball_r] = validate_parameters_L(outdir, data, dim_mfld, ball_r, params)
        
    validate_outdir(outdir);
    validate_data(data, dim_mfld);
        
    assert(isvector(ball_r) && all(isfinite(ball_r)) && isreal(ball_r) && all(ball_r>0) && length(ball_r)<=size(data,1));

    point_ids    = sort(params.Results.point_ids);
    [calib_points, reorder] = sort(params.Results.calib_points);    
    if (isempty(calib_points))        
        assert(isscalar(ball_r), 'When calib_points is unspecified, ball_r must be a scalar');        
    else
        assert(length(calib_points)==length(ball_r), 'Length of ball_r and calib_points must be equal');   
        ball_r = ball_r(reorder);
    end
    
end
