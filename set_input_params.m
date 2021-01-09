function params = set_input_params(data)
    
    params = inputParser;
    params.KeepUnmatched = false;
    params.CaseSensitive = true;
    params.PartialMatching = false;
     
    addParameter(params, 'calib_points',     [],                @(x) validate_point_list(x, size(data,1)));
    addParameter(params, 'point_ids',        [1:size(data,1)]', @(x) validate_point_list(x, size(data,1)));
    addParameter(params, 'global_r_prctile', NaN,               @(x) isscalar(x) && x>=0 && x<=100);

end
