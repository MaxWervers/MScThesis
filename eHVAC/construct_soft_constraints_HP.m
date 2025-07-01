function [E_in, e_in] = construct_soft_constraints_HP(N, mc, n, Tmin_air, Tmax_air, Tmax_w, u_prev_peak)
    % CONSTRUCT_SOFT_CONSTRAINTS Constructs inequality constraint matrix A and vector b
    % for the LP formulation with soft constraints on air temperature and peak demand.
    %
    % Inputs:
    %   N    - Number of time steps
    %   mc   - Number of control inputs (m)
    %   n    - Number of state variables
    %   Tmin - Nx1 vector of minimum temperature bounds
    %   Tmax - Nx1 vector of maximum temperature bounds
    %   u_prev_peak - Scalar value for previous peak demand
    %
    % Outputs:
    %   A - Constraint matrix ((3N) x (N*mc + N*n + N + 1))
    %   b - Constraint vector ((3N) x 1)

    % Index of T_air and T_water in state vector X
    ia = 1; % Adjust this if needed
    iw = 3; 

    % Compute sizes of different parts of the decision variable
    num_U = N * mc;  % Control input variables
    num_X = N * n;   % State variables
    num_zeta1 = N;   % Slack variable zeta1
    num_zeta2 = N;   % Slack variable zeta2
    num_z = 1;       % Peak demand variable

    % Construct selection matrix I_T1 and I_T3 to extract T_air and T_water from X
    I_T1 = zeros(N, num_X);
    for k = 1:N
        I_T1(k, (k-1)*n + ia) = 1; % Selecting T_air from X
    end
    
    I_T3 = zeros(N, num_X);
    for k = 1:N
        I_T3(k, (k-1)*n + iw) = 1; % Selecting T_air from X
    end

    % Construct identity matrix for slack variables
    I_zeta = eye(N);

    % Construct A matrix for the inequalities
    A_temp = [zeros(N, num_U), -I_T1,  -I_zeta, zeros(N,num_zeta2) , zeros(N,1);
         zeros(N, num_U),  I_T1, -I_zeta, zeros(N,num_zeta2), zeros(N,1);
         zeros(N, num_U), I_T3, zeros(N, num_zeta1), -I_zeta, zeros(N,1)];

    % Construct b vector
    b_temp = [-Tmin_air'; Tmax_air'; Tmax_w*ones(N,1)];

    
    % ---- Construct Peak Demand Constraint: u[k] - z ≤ u_prev_peak ----
    
    % Initialize peak demand constraint matrix
    A_peak = zeros(N, num_U + num_X + num_zeta1 + num_zeta2 + num_z);
    
    % Construct constraint: u[k] - z ≤ u_prev_peak
    for k = 1:N
        u_idx_start = (k - 1) * mc + 1;  % Start index for u[k]
        u_idx_end = k * mc;  % End index for u[k]
        
        % Assign coefficients for u[k]
        A_peak(k, u_idx_start:u_idx_end) = eye(mc); 
        
        % Assign coefficient for z (last column)
        A_peak(k, end) = -1;
    end
    
    % Right-hand side for peak constraint
    b_peak = u_prev_peak * ones(N, 1);

    % ---- Combine Constraints ----
    E_in = [A_temp; A_peak];
    e_in = [b_temp; b_peak];

end

