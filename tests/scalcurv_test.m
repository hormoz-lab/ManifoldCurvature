function tests = scalcurv_test
    tests = functiontests(localfunctions);
end

function test_D3(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    rng(917878);
    X = [rand(10000,3), rand(10000,2)];
    outpath = sprintf('%s/out/D3', folder);
    verifyWarningFree(testCase, @() manifold_curvature(outpath, X, 2, 0.5));    
    load(sprintf('%s/Curvature.mat', outpath));    
    verifyGreaterThan(testCase, sum((S-2*dS)<0 & (S+2*dS)>0)/sum(~isnan(S)), 0.95);
end

function test_S2(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    rng(72839);
    X = randn(10000, 3);
    X = X./sqrt(sum(X.^2,2));
    outpath = sprintf('%s/out/S2', folder);
    verifyWarningFree(testCase, @() manifold_curvature(outpath, X, 2, 0.017));    
    load(sprintf('%s/Curvature.mat', outpath));    
    verifyGreaterThan(testCase, sum((S-2*dS)<2 & (S+2*dS)>2)/sum(~isnan(S)), 0.7);
end

function test_glob_scale(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    L1 = load(sprintf('%s/out/S2/Curvature.mat', folder));
    sf = 2;
    outpath = sprintf('%s/out/S2_scale', folder);
    verifyWarningFree(testCase, @() manifold_curvature(outpath, sf*L1.data, L1.dim_mfld, L1.se_targ/sf));
    L2 = load(sprintf('%s/Curvature.mat', outpath));
    tol = 1e-15;
    verifyEqual(testCase, L1.N_neighbors, L2.N_neighbors);
    verifyLessThan(testCase, max(abs(L1.ball_r./L2.ball_r-1/sf)), tol);
    verifyLessThan(testCase, max(abs(L1.S     ./L2.S  - sf^2)), tol);
    verifyLessThan(testCase, max(abs(L1.dS    ./L2.dS - sf^2)), tol);
    verifyLessThan(testCase, max(abs(L1.gof   ./L2.gof- 1   )), tol);    
    verifyLessThan(testCase, max(cellfun(@(x,y) max(abs(x./y - sf  )), L1.hessmats, L2.hessmats)), tol);
    verifyLessThan(testCase, max(cellfun(@(x,y) max(abs(x./y - 1   )), L1.vars,     L2.vars)),     tol);
    verifyLessThan(testCase, max(cellfun(@(x,y) max(abs(x./y - sf^2)), L1.CovB,     L2.CovB)),     tol);
end
    
