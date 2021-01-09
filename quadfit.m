function [CovB, E, Sigma, hessmats, logL] = quadfit(dataT, dataN, pad_dim)

    if (nargin == 2)
        pad_dim = 0;
    end

    N = size(dataT,1);
    d_t = size(dataT,2);
    d_n = size(dataN,2);

    triuid = logical(triu(ones(d_t)));
    [r, c] = find(triuid);

    quadT = dataT(:,r).*dataT(:,c);   
    
    [beta,Sigma,E,CovB,logL] = mvregress([ones(N,1),quadT],dataN);
    
    CovB(1:size(beta,1):end,:) = [];                    % Remove constant term in CovB
    CovB(:,1:size(beta,1):end) = [];
    rep = 1+repmat(r~=c,d_n,1);                         % Remove factor of two in non-diagonal Hessian terms
    CovB = CovB./(rep.*rep');
    CovB = padarray(CovB, (size(beta,1)-1)*pad_dim*ones(1,2), 0, 'post');
    
    if (nargout >= 4)
        hessmats = zeros(size(beta,1)-1,d_n+pad_dim);       
        hessmats(:,1:d_n) = beta(2:end,:);                  % Ignore constant term
        hessmats(r~=c,:) = hessmats(r~=c,:)/2;              % Overestimated non-diagonal Hessian terms since we removed redundant x_ij in fit
    end

end