function N_min = get_min_points_for_regression(dim_mfld, dim_amb)

    % For regression need to compute coefficients for each x_ij, for 
    % 1 <= i <= j <= dim_mfld, so need at least nchoosek(dim_fld+1,2)
    % points. Also need to compute a coefficient on the constant term. 
    % Lastly, mvregress needs 1 point sample point than coefficients to
    % fit.

    N_min = max([nchoosek(2+dim_mfld-1,2), dim_amb])+2;
    
end