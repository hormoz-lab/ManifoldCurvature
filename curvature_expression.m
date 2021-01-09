function [scalF, Jf, Hf] = curvature_expression(d_t, d_n)

    assert(isscalar(d_t) && (d_t==floor(d_t)) && d_t >= 2);
    assert(isscalar(d_n) && (d_n==floor(d_n)) && d_n >= 1);

    hessmats = sym('h', [d_t d_t d_n]);
    for i = 1:d_n
        hvars{i} = nonzeros(triu(squeeze(hessmats(:,:,i))));
        hessmats(:,:,i) = tril(squeeze(hessmats(:,:,i)).',-1)+triu(squeeze(hessmats(:,:,i)));
    end
    
    h2 = 4*reshape(hessmats(:)*hessmats(:).', [size(hessmats), size(hessmats)]);
    if (d_n > 1)
        hc = sym(zeros(d_t,d_t,d_t,d_t));
        for i = 1:d_n
            hc = hc+squeeze(h2(:,:,i,:,:,i));
        end
    else
        hc = h2;
    end
    hc = simplify(hc);

    %permute the ids of hc to follow convention of R_ijkl
    ida = [1 3 2 4];
    idb = [3 2 1 4];
    hc_a = permute(hc,ida);
    hc_b = permute(hc,idb);

    riem = hc_a - hc_b; %Gauss formula for R_ijkl
    
    ricc = sym(zeros(d_t,d_t));
    for i = 1:d_t
        ricc = ricc+squeeze(riem(:,i,:,i));
    end
    
    scal = sym(zeros(1,1));
    for i = 1:d_t
        scal = scal+ricc(i,i);
    end

    scal  = simplify(scal);    
    scalF = matlabFunction(scal, 'Vars', {[hvars{:}]});
    
    if (nargout > 1)
        J  = simplify(gradient(scal, vertcat(hvars{:})));    
        Jf = matlabFunction(J, 'Vars', {[hvars{:}]});    
    end
    
    if (nargout > 2)
        H = simplify(hessian(scal, vertcat(hvars{:})));
        Hf = matlabFunction(H, 'Vars', {[hvars{:}]});
    end

        
end