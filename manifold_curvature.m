function manifold_curvature(outdir, data, dim_mfld, se_targ, varargin)
%manifold_curvature computes point-wise scalar curvatures on a dataset
%
%   manifold_curvature(OUTDIR, DATA, DIM_MFLD, SE_TARG) computes the 
%   scalar curvature for each point (row of DATA), assuming the manifold 
%   dimension is DIM_MFLD, and using a ball radius that achieves an
%   uncertainty of SE_TARG for the coefficients in the Second Fundamental
%   Form (2FF). The following files are saved to OUTDIR:
%
%   Curvature.mat - contains all variables including scalar curvatures,
%     standard errors, goodness-of-fits, Hessian matrices for the 2FF (upper 
%     triangular in column-major order), covariance matrices for 2FF
%     (upper triangular in column-major order), calibrated ball radii, 
%     and the Voronoi cells used to map calibration points to datapoints.
%
%   Calibration.txt - Running log of ball radius calibration step.
%
%   Log.txt - Running log of curvature computation step.
%
%   The code can optionally be invoked with some extra paramters:
%
%   manifold_curvature(..., 'point_ids', point_ids) only computes
%   curvatures on the subset of points specified by the vector point_ids. 
%   By default curvature is computed on all points.
%
%   manifold_curvature(..., 'calib_points', calib_points) only computes
%   ball radii for the set of points specified in vector calib_points.
%   Defaults to point_ids. All point_ids in the Voronoi cell of each
%   calib_point use the same ball radius for computing curvature.
%
%   manifold_curvature(..., 'global_r_prctile', p) computes curvatures 
%   using a constant ball size equal to the p-th percentile of ball radii
%   values computed over calib_points.
%
%   Examples:
%
%       % Compute curvature of 10K points drawn from S2 in R3
%       X = randn(10000,3);
%       X = X./vecnorm(X,2,2);
%       manifold_curvature('S2_run', X, 2, 0.01);
%
%       % Only compute curvature for every tenth point
%       manifold_curvature('S2_run', X, 2, 0.01, 'point_ids', [1:10:length(X)]');
%
%       % Only use first 100 points for ball radius calibration
%       manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]');
%
%       % Only use first 100 points for ball radius calibration, and use
%       % the median of these to compute curvatures for all data points
%       manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]', ...
%                          'global_r_prctile', 50);
% 
%   See also calibrate_neighborhoods, curvature_at_length_scale.
%
%   If you use this code, please cite:
%
%   D. Sritharan*, S. Wang*, S. Hormoz. "Computing the Riemannian 
%   curvature of image patch and single-cell RNA sequencing data 
%   manifolds using extrinsic differential geometry." bioRxiv(2021), 
%   https://doi.org/10.1101/2021.01.08.425885
%
%   Authors: Duluxan Sritharan & Shu Wang. Harvard Medical School.

  
    assert(nargin >= 4, 'Not all required input arguments supplied');
    params = set_input_params(data);    
    parse(params, varargin{:});
    [point_ids, calib_points, global_r_prctile] = validate_parameters(outdir, data, dim_mfld, se_targ, params);
    clear params;
    
    diary([tempname '.txt']);
    diary on;
    fprintf('Writing log to temporary location: %s\n', get(0,'DiaryFile'));    
    [ball_r, calib_points] = calibrate_neighborhoods(data, dim_mfld, se_targ, calib_points);
    diary off;
    if (~exist(outdir, 'dir'))        
        mkdir(outdir);
    end
    copyfile(get(0,'DiaryFile'), [outdir '/Calibration.txt']);    
    
    if (isnan(global_r_prctile))
        curvature_at_length_scale(outdir, data, dim_mfld, ball_r, calib_points, 'point_ids', point_ids)
    else
        curvature_at_length_scale(outdir, data, dim_mfld, prctile(ball_r, global_r_prctile), 'point_ids', point_ids);
    end
        
    save(sprintf('%s/Curvature.mat', outdir), 'se_targ', 'global_r_prctile', '-append');
    
end

