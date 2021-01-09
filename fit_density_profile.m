function [dens_power, rho, E, trial_dists, counts] = fit_density_profile(dist_to_center, N_min)

    assert(issorted(dist_to_center));
    
    min_dist = dist_to_center(CurvParams.init_redundancy_factor*N_min);    
    max_dist = dist_to_center(end);
    
    assert(min_dist > 0);
    
    N_scales = CurvParams.N_radius_grid;
    trial_dists = exp(linspace(log(min_dist),log(max_dist), N_scales))';       
    counts = arrayfun(@(i) sum(dist_to_center <= trial_dists(i)), [1:N_scales]');
    
    [dens_fit, ~, E] = regress(log(counts(2:N_scales-2)), [log(trial_dists(2:N_scales-2)) ones(N_scales-3,1)]);
    
    dens_power = dens_fit(1);
    rho = exp(dens_fit(2));
    
end