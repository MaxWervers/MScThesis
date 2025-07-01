function [N, N_vec_known, N_vec_ext, Steps, Ta, c, Tmin, Tmax] = sim_parameters_Nvec3(Nh_add, num_days, start, H_split, meta_data)
    % SIM_PARAMETERS - Define simulation parameters
    %
    %
    % Inputs:
    %    Nh - Prediction horizon, in hours
    %    num_days - Number of days in simulation
    %    start - Start index related to start date
    %    H_split - Split hour into intervals
    %    meta_data - Metadata containing ambient temperatures, electricity prices, minimum and maximum temperatures
    %
    % Outputs:
    %    N - Prediction horizon, in steps
    %    Steps - Total simulation steps
    %    Ta - Ambient temperatures
    %    c - Electricity prices
    %    Tmin - Minimum temperatures
    %    Tmax - Maximum temperatures
    
    Nh_max = 35 + Nh_add; % Max amount of information on ce known

    % Compute total simulation steps
    H_sim = num_days * 24; % Total hours in simulation
    Steps = H_sim * H_split; % Total simulation steps
    
    % Compute prediction horizon in steps
    N = Nh_max * H_split;
    
    % Obtain and expand ambient temperatures
    Ta = meta_data{start:start+H_sim , 2}';
    Ta = repelem(Ta, H_split);
    
    % Obtain and expand electricity prices
    c = meta_data{start:start+H_sim , 3}';
    c = repelem(c, H_split);
    
    % Obtain and expand minimum temperatures
    Tmin = meta_data{start:start+H_sim , 4}';
    Tmin = repelem(Tmin, H_split);
    
    % Obtain and expand maximum temperatures
    Tmax = meta_data{start:start+H_sim , 5}';
    Tmax = repelem(Tmax, H_split);

    % Obtain prediction horizon and expand to correct split
    N_vec_small = meta_data{start:start+H_sim-1 , 6}' * H_split;

    N_vec_full = [];
    N_inter_add = linspace(0,H_split-1,H_split);

    for i = 1:length(N_vec_small)
        N_vec_full_loop = N_vec_small(i) - N_inter_add;
        N_vec_full = [N_vec_full, N_vec_full_loop];
    end
    
    N_vec_known = N_vec_full;
    N_vec_ext = N_vec_known + Nh_add*H_split;


    % Define how many final steps should shrink
    K_shrink = N_vec_ext(1);
    N_min = 1;     % Minimum horizon at final step (e.g., 1 hour)

    if Steps > K_shrink
    % Indices where the shrinking starts
    idx_shrink_start = Steps - K_shrink ;
    
    % Create shrinking horizon vector (linearly from current to N_min)
    shrink_horizon = linspace(N_vec_ext(idx_shrink_start), N_min, K_shrink+1);
    
    % Replace the tail of N_vec_ext
    N_vec_ext(idx_shrink_start:end) = shrink_horizon;
    end

    % Define how many final steps should shrink
    K_shrink2 = N_vec_known(1);
    N_min = 1;     % Minimum horizon at final step (e.g., 1 hour)

    if Steps > K_shrink
    % Indices where the shrinking starts
    idx_shrink_start2 = Steps - K_shrink2 ;
    
    % Create shrinking horizon vector (linearly from current to N_min)
    shrink_horizon2 = linspace(N_vec_known(idx_shrink_start2), N_min, K_shrink2+1);
    
    % Replace the tail of N_vec_ext
    N_vec_known(idx_shrink_start2:end) = shrink_horizon2;
    end
end

