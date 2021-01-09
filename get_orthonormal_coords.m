function [local_coords, dim_nt, vars] = get_orthonormal_coords(data)

    [~, local_coords, ~, ~, vars] = pca(data);
    if (~isempty(vars))
        dim_nt = find(vars/vars(1) > CurvParams.trivial_dim_var, 1, 'last');
    else
        dim_nt = 0;
    end
    
end