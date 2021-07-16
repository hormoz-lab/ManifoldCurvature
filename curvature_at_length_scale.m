function curvature_at_length_scale(outdir, data, dim_mfld, ball_r, varargin)
%curvature_at_length_scale computes point-wise scalar curvatures on a dataset
%at the specified length scale
%
%   curvature_at_length_scale(OUTDIR, DATA, DIM_MFLD, BALL_R) computes the 
%   scalar curvature for each point (row of DATA), assuming the manifold 
%   dimension is DIM_MFLD, and using a constant ball radius given by the 
%   scalar BALL_R.
%
%   curvature_at_length_scale(OUTDIR, DATA, DIM_MFLD, BALL_R, CALIB_POINTS) 
%   allows variable ball radii given by a vector BALL_R corresponding to
%   CALIB_POINTS. All datapoints (rows of DATA) in the Voronoi cell of each
%   CALIB_POINT will use the same BALL_R.
%
%   The following files are saved to OUTDIR:
%
%   Curvature.mat - contains all variables including scalar curvatures,
%     standard errors, goodness-of-fits, Hessian matrices for the 2FF (upper 
%     triangular in column-major order), covariance matrices for 2FF
%     (upper triangular in column-major order), calibrated ball radii, 
%     and the Voronoi cells used to map calibration points to datapoints.
%
%   Log.txt - Running log of curvature computation step.
%
%   The code can optionally be invoked with some extra paramters:
%
%   curvature_at_length_scale(..., 'point_ids', point_ids) only computes
%   curvatures on the subset of points specified by the vector point_ids. 
%   By default curvature is computed on all points.
%
%   Examples:
%
%       % Curvature of 10K points on S2 in R3 using constant ball radius of 0.2
%       X = randn(10000,3);
%       X = X./vecnorm(X,2,2);
%       curvature_at_length_scale('S2_run', X, 2, 0.2);
%
%       % Only compute curvature for every tenth point
%       curvature_at_length_scale('S2_run', X, 2, 0.2, 'point_ids', [1:10:length(X)]');
%
%       % Use variable ball radii sizes
%       [R, ids] = calibrate_neighborhoods(X, 2, 0.01, [1:100]');
%       curvature_at_length_scale('S2_run', X, 2, R, ids);
%
%       % This is equivalent to
%       manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]');
% 
%   See also manifold_curvature, calibrate_neighborhoods.
%
%   If you use this code, please cite:
%
%   D. Sritharan*, S. Wang*, S. Hormoz. "Computing the Riemannian 
%   curvature of image patch and single-cell RNA sequencing data 
%   manifolds using extrinsic differential geometry." PNAS (2021),
%   https://doi.org/10.1073/pnas.2100473118
%
%   Authors: Duluxan Sritharan & Shu Wang. Harvard Medical School.


    assert(nargin >= 4, 'Not all required input arguments supplied');
    params = set_input_params_L(data);
    parse(params, varargin{:});
    [point_ids, calib_points, ball_r] = validate_parameters_L(outdir, data, dim_mfld, ball_r, params);
    clear params;
 
    tic;    
    diary([tempname '.txt']);
    diary on;
    fprintf('Writing log to temporary location: %s\n', get(0,'DiaryFile'));
    
    dim_amb = size(data,2);
    N_min = get_min_points_for_regression(dim_mfld, dim_amb);    
    
    fprintf('Generating curvature expression...\n');
    [scalF, Jf] = curvature_expression(dim_mfld, dim_amb-dim_mfld);
   
    fprintf('Computing curvature...\n');
    
    N           = length(point_ids);
    S           = NaN(N,1);    
    dS          = NaN(N,1);
    gof         = NaN(N,1);
    N_neighbors = NaN(N,1);
    hessmats    = cell(N,1);
    vars        = cell(N,1);
    CovB        = cell(N,1);    
   
    if (isscalar(ball_r))
        which_calib = ones(N,1);
    else
        which_calib = knnsearch(data(calib_points,:), data(point_ids,:));        
    end
    
    spmd            
        warning('off', 'stats:pca:ColRankDefX');
        warning('off', 'MATLAB:nearlySingularMatrix');
        warning('off', 'MATLAB:rankDeficientMatrix');        
    end
       
    parfor j = 1:N
        [~, neighbor_id] = pdist2(data, data(point_ids(j),:), 'Euclidean', 'Radius', ball_r(which_calib(j))+CurvParams.pdist_eps);        
        neighbor_id = neighbor_id{1};
        N_neighbors(j) = length(neighbor_id);        
        if (N_neighbors(j) >= CurvParams.fit_redundancy_factor*N_min)            
            try
                [CovB{j}, gof(j), S(j), dS(j), hessmats{j}, vars{j}] = neighborhood_compute(data(neighbor_id,:), dim_mfld, scalF, Jf);                         
                fprintf('Computed point %10d/%d: S=%10.5f dS=%10.5f gof=%10.9f r=%10.5f N=%10d\n', ...
                    j, N, S(j), dS(j), gof(j), ball_r(which_calib(j)), N_neighbors(j));
            catch e
                fprintf('Computed point %10d/%d: FAILED - %s: %s\n', j, N, e.identifier, e.message);
            end            
        end
    end
    
    fprintf('...mean+/SD number of neighbours = (%g,%g)\n', nanmean(N_neighbors), nanstd(N_neighbors));
    fprintf('...successfully computed curvature for %d/%d points\n', sum(~isnan(S)), N);
       
    if (~exist(outdir, 'dir'))        
        mkdir(outdir);
    end
    
    try    
        save(sprintf('%s/Curvature.mat', outdir), 'data', 'dim_mfld', 'point_ids', 'ball_r', 'calib_points', 'which_calib', ...
            'S', 'dS', 'gof', 'N_neighbors', 'hessmats', 'CovB', 'vars');
    catch
        save(sprintf('%s/Curvature.mat', outdir), 'data', 'dim_mfld', 'point_ids', 'ball_r', 'calib_points', 'which_calib', ...
            'S', 'dS', 'gof', 'N_neighbors', 'hessmats', 'CovB', 'vars', '-v7.3');            
    end
        
    fprintf('Completed in %6.4g seconds\n', toc);

    diary off;
    copyfile(get(0,'DiaryFile'), [outdir '/Log.txt']);
    
end

