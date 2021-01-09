function validate_point_list(x, N)

    assert(isvector(x) && all(isfinite(x)) && length(unique(x))==length(x) && all(x>0 & x<=N));
    
end
