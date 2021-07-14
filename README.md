# Manifold Curvature

This code computes the point-wise scalar curvature for a set of datapoints. It was written in MATLAB 2019a and may not be compatible with earlier versions.

## Citation

Please cite the [companion paper](https://doi.org/10.1101/2021.01.08.425885) if you use this software package:

> D. Sritharan*, S. Wang*, S. Hormoz. "Computing the Riemannian curvature of image patch and single-cell RNA sequencing data manifolds using extrinsic differential geometry." bioRxiv(2021).

## Download

To download the code into the directory given by CODE_PATH:

```bash
$ git clone https://gitlab.com/hormozlab/ManifoldCurvature.git CODE_PATH
```

## Installation

Add the ManifoldCurvature software to the MATLAB path:
	
```MATLAB
>> addpath(genpath(CODE_PATH), '-end');
```

If you are installing on your local machine, you should additionally save the path.

```MATLAB
>> savepath;
```

If you are installing on a cluster machine, the MATLAB path file will be read-only, so you will have to re-run the first command at the start of each session.

## Usage

The main function is `manifold_curvature`, which computes scalar curvatures at a different length scale for each point so that the uncertainty in the Second Fundamental Form coefficients are upper-bounded by `se_targ`.

```MATLAB
>> help manifold_curvature
```

    manifold_curvature computes point-wise scalar curvatures on a dataset
 
    manifold_curvature(OUTDIR, DATA, DIM_MFLD, SE_TARG) computes the 
    scalar curvature for each point (row of DATA), assuming the manifold 
    dimension is DIM_MFLD, and using a ball radius that achieves an
    uncertainty of SE_TARG for the coefficients in the Second Fundamental
    Form (2FF). The following files are saved to OUTDIR:
 
    Curvature.mat - contains all variables including scalar curvatures,
      errorbars, goodness-of-fits, Hessian matrices for the 2FF (upper 
      triangular in column-major order), covariance matrices for 2FF
      (upper triangular in column-major order), calibrated ball radii, 
      and the Voronoi cells used to map calibration points to datapoints.
 
    Calibration.txt - Running log of ball radius calibration step.
 
    Log.txt - Running log of curvature computation step.
 
    The code can optionally be invoked with some extra paramters:
 
    manifold_curvature(..., 'point_ids', point_ids) only computes
    curvatures on the subset of points specified by the vector point_ids. 
    By default curvature is computed on all points.
 
    manifold_curvature(..., 'calib_points', calib_points) only computes
    ball radii for the set of points specified in vector calib_points.
    Defaults to point_ids. All point_ids in the Voronoi cell of each
    calib_point use the same ball radius for computing curvature.
 
    manifold_curvature(..., 'global_r_prctile', p) computes curvatures 
    using a constant ball size equal to the p-th percentile of ball radii
    values computed over calib_points.
 
    Examples:
 
        % Compute curvature of 10K points drawn from S2 in R3
        X = randn(10000,3);
        X = X./vecnorm(X,2,2);
        manifold_curvature('S2_run', X, 2, 0.01);
 
        % Only compute curvature for every tenth point
        manifold_curvature('S2_run', X, 2, 0.01, 'point_ids', [1:10:length(X)]');
 
        % Only use first 100 points for ball radius calibration
        manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]');
 
        % Only use first 100 points for ball radius calibration, and use
        % the median of these to compute curvatures for all data points
        manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]', ...
                           'global_r_prctile', 50);
  
    See also calibrate_neighborhoods, curvature_at_length_scale.
 
    Authors: Duluxan Sritharan & Shu Wang. Harvard Medical School.

To choose `se_targ`, use the `calibrate_neighborhoods` function and verify if the distribution of GOF p-values is flat and the neighborhood sizes are appropriate given the global length scale of the data.

```MATLAB
>> help calibrate_neighborhoods
```

    calibrate_neighborhoods computes point-wise ball radii for regression
 
    [BALL_R, CALIB_POINTS, GOF, N_NEIGHBORS] = ...
    calibrate_neighborhoods(DATA, DIM_MFLD, SE_TARG) computes the ball 
    radius at each point (row of DATA) needed to achieve a target 
    uncertainty of SE_TARG for the diagonal coefficients in the Second 
    Fundamental Form, assuming the manifold dimension is DIM_MFLD. 
    CALIB_POINTS is a vector of point IDs for which calibration was 
    successful, and BALL_R is a vector of corresponding ball radii. The
    function also returns goodness-of-fits, GOF, and the number of points 
    in the ball, N_NEIGHBORS, for each calibration point.
 
    [BALL_R, CALIB_POINTS, GOF, N_NEIGHBORS] = ...
    calibrate_neighborhoods(DATA, DIM_MFLD, SE_TARG, CALIB_POINTS) 
    only calibrates for the set of points specified in vector CALIB_POINTS.
 
    Examples:
 
        % Compute ball radii for 10K points drawn from S2 in R3
        X = randn(10000,3);
        X = X./vecnorm(X,2,2);
        [R, ids] = calibrate_neighborhoods(X, 2, 0.01);
 
        % Only calibrate first 100 points;
        [R, ids] = calibrate_neighborhoods(X, 2, 0.01, [1:100]');
 
    See also manifold_curvature, calibrate_neighborhoods.
 
    Authors: Duluxan Sritharan & Shu Wang. Harvard Medical School.

To compute curvatures at a specific length scale for each point use the `curvature_at_length_scale` function.

```MATLAB
>> help curvature_at_length_scale
```

    curvature_at_length_scale computes point-wise scalar curvatures on a dataset
    at the specified length scale
 
    curvature_at_length_scale(OUTDIR, DATA, DIM_MFLD, BALL_R) computes the 
    scalar curvature for each point (row of DATA), assuming the manifold 
    dimension is DIM_MFLD, and using a constant ball radius given by the 
    scalar BALL_R.
 
    curvature_at_length_scale(OUTDIR, DATA, DIM_MFLD, BALL_R, CALIB_POINTS) 
    allows variable ball radii given by a vector BALL_R corresponding to
    CALIB_POINTS. All datapoints (rows of DATA) in the Voronoi cell of each
    CALIB_POINT will use the same BALL_R.
 
    The following files are saved to OUTDIR:
 
    Curvature.mat - contains all variables including scalar curvatures,
      errorbars, goodness-of-fits, Hessian matrices for the 2FF (upper 
      triangular in column-major order), covariance matrices for 2FF
      (upper triangular in column-major order), calibrated ball radii, 
      and the Voronoi cells used to map calibration points to datapoints.
 
    Log.txt - Running log of curvature computation step.
 
    The code can optionally be invoked with some extra paramters:
 
    curvature_at_length_scale(..., 'point_ids', point_ids) only computes
    curvatures on the subset of points specified by the vector point_ids. 
    By default curvature is computed on all points.
 
    Examples:
 
        % Curvature of 10K points on S2 in R3 using constant ball radius of 0.2
        X = randn(10000,3);
        X = X./vecnorm(X,2,2);
        curvature_at_length_scale('S2_run', X, 2, 0.2);
 
        % Only compute curvature for every tenth point
        curvature_at_length_scale('S2_run', X, 2, 0.2, 'point_ids', [1:10:length(X)]');
 
        % Use variable ball radii sizes
        [R, ids] = calibrate_neighborhoods(X, 2, 0.01, [1:100]');
        curvature_at_length_scale('S2_run', X, 2, R, ids);
 
        % This is equivalent to
        manifold_curvature('S2_run', X, 2, 0.01, 'calib_points', [1:100]');
  
    See also manifold_curvature, calibrate_neighborhoods.

    Authors: Duluxan Sritharan & Shu Wang. Harvard Medical School.

## Troubleshooting

If you're running into issues, check if the ManifoldCurvature code runs properly on test datasets on your installation.

```MATLAB
>> cd([CODE_PATH '/tests']);
>> runtests;
```

#### Prepared By: Duluxan Sritharan
