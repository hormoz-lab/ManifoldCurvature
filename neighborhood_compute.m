function [CovB, gof, S, dS, hessmats, vars] = neighborhood_compute(data, dim_mfld, scalF, Jf)

    if (nargin == 2)
        calibrate_mode = true;
        assert(nargout <= 2);
    elseif (nargin == 4)
        calibrate_mode = false;
    else
        error('Incorrect argument list');
    end  

    dim_amb = size(data,2);
    [local_coords, dim_nt, vars] = get_orthonormal_coords(data);        
    
    % Try to fit at least one normal coordinate
    dim_nt = max(dim_nt,dim_mfld+1);
    if (calibrate_mode == false)
        [CovB, E, Sigma, hessmats] = quadfit(local_coords(:,1:dim_mfld), local_coords(:,(dim_mfld+1):dim_nt), dim_amb-dim_nt);
        h = hessmats(:);
        S = scalF(h);
        Jfval = Jf(h);
        dS = sqrt(Jfval'*CovB*Jfval);
        CovB = CovB(triu(true(size(CovB))));
    else
        [CovB, E, Sigma]           = quadfit(local_coords(:,1:dim_mfld), local_coords(:,(dim_mfld+1):dim_nt), dim_amb-dim_nt);
    end
    if (nargout >= 2)
        gof = compute_gof(E, Sigma);
    end
  
end
                
                
                
                
                
                