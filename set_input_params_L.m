function params = set_input_params_L(data)
    
    params = inputParser;
    params.KeepUnmatched = false;
    params.CaseSensitive = true;
    params.PartialMatching = false;
   
    addOptional (params, 'calib_points', [],                @(x) validate_point_list(x, size(data,1)));
    addParameter(params, 'point_ids',    [1:size(data,1)]', @(x) validate_point_list(x, size(data,1)));        
    
end
