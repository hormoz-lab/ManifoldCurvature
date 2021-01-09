function [ball_r, se, gof, N_ball, num_iter, dens_power, rho] = calibrate_at_point(data, calib_point, dim_mfld, se_targ)
    
    dim_amb = size(data,2);
    N_min = get_min_points_for_regression(dim_mfld, dim_amb);    
    [dist_to_center, point_id] = sort(pdist2(data,data(calib_point,:)));    
    [dens_power, rho] = fit_density_profile(dist_to_center, N_min);    
    
    warning('off', 'stats:pca:ColRankDefX');
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    
    N_ball_list = zeros(CurvParams.max_iter,1);
    num_iter = 0;
    ball_r = dist_to_center(ceil(CurvParams.init_redundancy_factor*N_min));
    
    while (num_iter<CurvParams.max_iter && ball_r < dist_to_center(end))
        
        num_iter = num_iter+1;
        
        N_ball = sum(dist_to_center <= ball_r);
        if ((N_ball<N_min) || ismember(N_ball, N_ball_list))
            break;
        end
        
        CovB = neighborhood_compute(data(point_id(1:N_ball),:), dim_mfld);
        se = sqrt(max(abs(diag(CovB))));
        alpha = se/se_targ;
        
        if (abs(alpha-1)<CurvParams.se_targ_tol)
            break;
        else
            ball_r = ball_r*(alpha)^(2/(dens_power+4));
            % Err on the bigger side if ball size has to be increased
            if (alpha>1)
                if (ball_r < dist_to_center(N_ball+1))
                    ball_r = dist_to_center(N_ball+1);                   
                end
            else                
                if (ball_r > dist_to_center(N_ball-1))                    
                    break;                    
                end
            end
            N_ball_list(num_iter) = N_ball;
        end       
    end
    
    ball_r = min(max(ball_r, dist_to_center(ceil(CurvParams.fit_redundancy_factor*N_min))), dist_to_center(end));
    N_ball = sum(dist_to_center <= ball_r);
    [CovB, gof] = neighborhood_compute(data(point_id(1:N_ball),:), dim_mfld);    
    se = sqrt(max(abs(diag(CovB))));    
    
end
    
    