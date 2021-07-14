function validate_data(data, dim_mfld)

    assert(all(isfinite(data(:))), 'data should be a matrix of finite scalar values');
    assert(isreal(data), 'data should be a matrix of real values');
    assert(size(data,2)>2, 'Ambient dimension needs to be > 2');
    assert(isscalar(dim_mfld) && isreal(dim_mfld) && dim_mfld<size(data,2) && dim_mfld>=2, 'dim_mfld should be scalar in [2, size(data,2)]');
    
    [~, ~, which] = unique(data, 'rows');
    ind = arrayfun(@(i) mat2str(find(which==i)), find(accumarray(which,1)>1), 'un', false);
    assert(isempty(ind), 'There should be no duplicate points. The following point_ids are duplicates: %s\n', horzcat(ind{:}))
    
    N_data_min = get_min_points_for_regression(dim_mfld, size(data,2))*CurvParams.init_redundancy_factor;
    assert(size(data,1)>= N_data_min, sprintf('Need at least %d data points for manifold dimension %d', N_data_min, dim_mfld));
    
end