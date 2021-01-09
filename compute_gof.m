function gof = compute_gof(E, Sigma)
    
    [N, n] = size(E);    
    E = (E-mean(E,1));
    temp = sum(sum((E*inv(Sigma)).*E,2).^2)/N;
    temp = sqrt(N/(8*n*(n+2)))*(temp-n*(n+2));
    gof = 2*(1-normcdf(abs(temp)));
    
end