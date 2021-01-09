function tests = pipeline_test
    tests = functiontests(localfunctions);
end

function test_dynamic_neighborhoods(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    rng(49954);
    X = gamrnd(10, 5, [1000, 7]);
    outpath = sprintf('%s/out/R_dynam', folder);
    verifyWarningFree(testCase, @() manifold_curvature(outpath, X, 4, 0.05, 'calib_points', [1:100]', 'point_ids', [1:2:1000]'));
    check_all_files(testCase, outpath, true);
end

function test_specified_neighborhoods(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));  
    rng(49954);
    X = gamrnd(10, 5, [1000, 7]);    
    output_path = sprintf('%s/out/R_fixed', folder);
    [R, ids] = verifyWarningFree(testCase, @() calibrate_neighborhoods(X, 4, 0.05, [1:100]'));
    verifyWarningFree(testCase, @() curvature_at_length_scale(output_path, X, 4, R, ids, 'point_ids', [1:2:1000]'));
    check_all_files(testCase, output_path, false);
end

function test_equiv_pipeline(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath')); 
    dynam = load(sprintf('%s/out/R_dynam/Curvature.mat', folder));
    fixed = load(sprintf('%s/out/R_fixed/Curvature.mat',  folder));
    verifyEqual(testCase, dynam.data,         fixed.data        );
    verifyEqual(testCase, dynam.dim_mfld,     fixed.dim_mfld    );
    verifyEqual(testCase, dynam.ball_r,       fixed.ball_r      );
    verifyEqual(testCase, dynam.point_ids,    fixed.point_ids   );
    verifyEqual(testCase, dynam.calib_points, fixed.calib_points);
    verifyEqual(testCase, dynam.S,            fixed.S           );
    verifyEqual(testCase, dynam.dS,           fixed.dS          );
    verifyEqual(testCase, dynam.gof,          fixed.gof         );
    verifyEqual(testCase, dynam.N_neighbors,  fixed.N_neighbors );
    verifyEqual(testCase, dynam.hessmats,     fixed.hessmats    );
    verifyEqual(testCase, dynam.vars,         fixed.vars        );
    verifyEqual(testCase, dynam.CovB,         fixed.CovB        );
    verifyEqual(testCase, dynam.which_calib,  fixed.which_calib );
    verifyEqual(testCase, dynam.global_r_prctile, NaN);
end

function test_const_neighborhoods(testCase)
    [folder, ~, ~] = fileparts(mfilename('fullpath'));
    rng(37754);
    X = poissrnd(12, [1000, 8]);
    outpath = sprintf('%s/out/R_const', folder);
    verifyWarningFree(testCase, @() manifold_curvature(outpath, X, 3, 0.05, 'calib_points', [1:100]', 'point_ids', [1:2:1000]', 'global_r_prctile', 50));
    check_all_files(testCase, outpath, true);
    [folder, ~, ~] = fileparts(mfilename('fullpath')); 
    const = load(sprintf('%s/out/R_const/Curvature.mat', folder));
    verifyTrue(testCase, all(const.which_calib==1));
    verifyTrue(testCase, isscalar(const.ball_r));
    verifyFalse(testCase, isnan(const.global_r_prctile));
end

function check_all_files(testCase, outpath, full)
    verifyEqual(testCase, exist(sprintf('%s/Curvature.mat',       outpath), 'file'), 2);    
    verifyEqual(testCase, exist(sprintf('%s/Log.txt',             outpath), 'file'), 2);
    load(sprintf('%s/Curvature.mat', outpath));
    n = size(data,2);
    N = length(point_ids);    
    c = max(length(calib_points), 1);
    verifyEqual(testCase, size(ball_r), [c, 1]);
    verifyEqual(testCase, size(S), [N, 1]);
    verifyEqual(testCase, size(dS), [N, 1]);
    verifyEqual(testCase, size(gof), [N, 1]);
    verifyEqual(testCase, size(N_neighbors), [N, 1]);
    verifyEqual(testCase, size(which_calib), [N, 1]);    
    verifyEqual(testCase, size(hessmats), [N, 1]);
    verifyEqual(testCase, size(vars), [N, 1]);
    verifyEqual(testCase, size(CovB), [N, 1]);    
    mask = ~isnan(S);
    verifyTrue(testCase, all(cellfun(@(x) isequal(size(x), [n, 1]), vars(mask))));
    codim = n-dim_mfld;
    coeffs_per_codim = dim_mfld*(dim_mfld+1)/2;
    verifyTrue(testCase, all(cellfun(@(x) isequal(size(x), [coeffs_per_codim, codim]), hessmats(mask))));
    uniq_coeffs = codim * coeffs_per_codim;
    verifyTrue(testCase, all(cellfun(@(x) isequal(size(x), [uniq_coeffs*(uniq_coeffs+1)/2, 1]), CovB(mask))));
    if (full)
        verifyEqual(testCase, exist(sprintf('%s/Calibration.txt', outpath), 'file'), 2);
        verifyTrue(testCase, isscalar(se_targ));
        verifyTrue(testCase, isscalar(global_r_prctile));
    else    
        verifyFalse(testCase, exist(sprintf('%s/Calibration.txt', outpath), 'file')==2);
        verifyFalse(testCase, exist('se_targ')==1);
        verifyFalse(testCase, exist('global_r_prctile')==1);
    end
end
