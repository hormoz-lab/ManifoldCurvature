classdef (Sealed) CurvParams < handle
        
    %% Parameters for program - put constants/globals here
    properties (Constant, GetAccess=public)
        
        init_redundancy_factor = 2.0;   % N_min times this many points for initial guess for ball radius
        fit_redundancy_factor = 2.0;    % N_min times this many points needed when computing curvature        
        N_radius_grid = 10;             % Number of grid points over which to perform density fitting
        ball_r_tol = 1e-6;              % Convergence tolerance for r and r_guess
        se_targ_tol = 1e-2;             % Fractional convergence tolerance for se_targ        
        pdist_eps = 1e-15;              % Margin since pdist computes < r and we want <=r (machine epsilon for double-precision)
        trivial_dim_var = 1e-9;         % Relative variance of trivial dimension
        max_iter = 10;                  % Maximum iterations for neighborhood ball size
    end
   
end
