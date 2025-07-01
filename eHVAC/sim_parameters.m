function [N, Steps, Ta, c, Tmin, Tmax] = sim_parameters(Nh, num_days, start, H_split, meta_data)
    % SIM_PARAMETERS - Define simulation parameters
    %
    % Syntax:  [N, Steps, Ta, c, Tmin, Tmax] = sim_parameters(Nh, num_days, start, H_split, meta_data)
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
    
    % Compute total simulation steps
    H_sim = num_days * 24; % Total hours in simulation
    Steps = H_sim * H_split; % Total simulation steps
    
    % Compute prediction horizon in steps
    N = Nh * H_split;
    
    % Obtain and expand ambient temperatures
    Ta = meta_data{start:start+H_sim+Nh-1 , 2}';
    Ta = repelem(Ta, H_split);
    
    % Obtain and expand electricity prices
    c = meta_data{start:start+H_sim+Nh-1 , 3}';
    c = repelem(c, H_split);
    
    % Obtain and expand minimum temperatures
    Tmin = meta_data{start:start+H_sim+Nh-1 , 4}';
    Tmin = repelem(Tmin, H_split);
    
    % Obtain and expand maximum temperatures
    Tmax = meta_data{start:start+H_sim+Nh-1 , 5}';
    Tmax = repelem(Tmax, H_split);
end

