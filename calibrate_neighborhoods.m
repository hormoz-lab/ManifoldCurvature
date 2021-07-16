function [ball_r, calib_points, gof, N_neighbors] = calibrate_neighborhoods(data, dim_mfld, se_targ, calib_points) 
%calibrate_neighborhoods computes point-wise ball radii for regression
%
%   [BALL_R, CALIB_POINTS, GOF, N_NEIGHBORS] = ...
%   calibrate_neighborhoods(DATA, DIM_MFLD, SE_TARG) computes the ball 
%   radius at each point (row of DATA) needed to achieve a target 
%   uncertainty of SE_TARG for the diagonal coefficients in the Second 
%   Fundamental Form, assuming the manifold dimension is DIM_MFLD. 
%   CALIB_POINTS is a vector of point IDs for which calibration was 
%   successful, and BALL_R is a vector of corresponding ball radii. The
%   function also returns goodness-of-fits, GOF, and the number of points 
%   in the ball, N_NEIGHBORS, for each calibration point.
%
%   [BALL_R, CALIB_POINTS, GOF, N_NEIGHBORS] = ...
%   calibrate_neighborhoods(DATA, DIM_MFLD, SE_TARG, CALIB_POINTS) 
%   only calibrates for the set of points specified in vector CALIB_POINTS.
%
%   Examples:
%
%       % Compute ball radii for 10K points drawn from S2 in R3
%       X = randn(10000,3);
%       X = X./vecnorm(X,2,2);
%       [R, ids] = calibrate_neighborhoods(X, 2, 0.01);
%
%       % Only calibrate first 100 points;
%       [R, ids] = calibrate_neighborhoods(X, 2, 0.01, [1:100]');
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
    

    assert(nargin >=3, 'Not all required input arguments supplied');
    validate_data(data, dim_mfld);
    validate_se(se_targ);
    
    if (nargin == 3)
        calib_points = [1:size(data,1)]';
    else
        validate_point_list(calib_points, size(data,1));
    end
    
    tic;
    fprintf('Calibrating neighborhood sizes...\n');
    
    N_calib = length(calib_points);    
    
    ball_r          = NaN(N_calib,1);
    se_achieved     = NaN(N_calib,1);
    gof             = NaN(N_calib,1);
    N_neighbors     = NaN(N_calib,1);
    
    parfor t = 1:N_calib        
        try
            [ball_r(t), se_achieved(t), gof(t), N_neighbors(t), num_iter] = calibrate_at_point(data, calib_points(t), dim_mfld, se_targ);
            fprintf('Calibrated point %10d/%d: r=%10.5f se=%10.5f gof=%10.9f N_ball=%10d iters=%3d\n', ...
                t, N_calib, ball_r(t), se_achieved(t), gof(t), N_neighbors(t), num_iter);
        catch e
            fprintf('Calibrated point %10d/%d: FAILED - %s: %s\n', t, N_calib, e.identifier, e.message);
        end
    end    
    
    mask = ~isnan(ball_r);
    calib_points = calib_points(mask);
    ball_r       = ball_r      (mask);
    se_achieved  = se_achieved (mask);
    gof          = gof         (mask);
    N_neighbors  = N_neighbors (mask);
    
    fprintf('...Target SE: %6.4g, Achieved SE: %6.4g +- %6.4g\n', se_targ, median(se_achieved), std(se_achieved));    
    fprintf('...Average ball radius = %10.5g\n', mean(ball_r));
    fprintf('...Calibration successful for %d/%d points\n', length(ball_r), N_calib);
    
    fprintf('Completed in %6.4g seconds\n', toc);
    
end