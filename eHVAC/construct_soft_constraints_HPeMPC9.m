function [E_in, e_in] = construct_soft_constraints_HPeMPC9(N, mc, n, Tmin_air, Tmax_air, Tmax_w, u_prev_peak, z_prev, a_i, b_i)
    % CONSTRUCT_SOFT_CONSTRAINTS_HPEMPC6 - Construct inequality constraints for HPeMPC6
    %
    % Inputs:
    %    N - Prediction horizon, in steps (scalar)
    %    mc - Number of control inputs
    %    n - Number of states
    %    Tmin_air - Lower bound for air temperature (N x 1)
    %    Tmax_air - Upper bound for air temperature (N x 1)
    %    Tmax_w - Upper bound for water temperature (scalar)
    %    u_prev_peak - Previous peak demand (scalar)
    %    a_i - Coefficients for terminal cost constraints (num_i x 1)
    %    b_i - Right-hand side for terminal cost constraints (num_i x 1)
    %
    % Outputs:
    %    E_in - Inequality matrix for the optimization problem
    %    e_in - Right-hand side for the inequality constraints

    % Index of T_air and T_water in state vector X
    ia = 1; 
    iw = 3; 

    % Compute sizes of different parts of the decision variable
    num_U = N * mc;  % Control input variables
    num_X = N * n;   % State variables
    num_zeta1 = N;   % Slack variable zeta1
    num_zeta2 = N;   % Slack variable zeta2
    num_i = length(b_i); % Amount off constraints due to terminal cost
    num_z = 1;       % Peak demand variable
    num_v = 1;

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
    A_temp = [zeros(N, num_U), -I_T1,  -I_zeta, zeros(N,num_zeta2) , zeros(N,1), zeros(N,1);
         zeros(N, num_U),  I_T1, -I_zeta, zeros(N,num_zeta2), zeros(N,1), zeros(N,1);
         zeros(N, num_U), I_T3, zeros(N, num_zeta1), -I_zeta, zeros(N,1), zeros(N,1)];

    % Construct b vector
    b_temp = [-Tmin_air'; Tmax_air'; Tmax_w*ones(N,1)];

    
    % ---- Construct Peak Demand Constraint: u[k] - z ≤ u_prev_peak ----
    
    % Initialize peak demand constraint matrix
    A_peak = zeros(N, num_U + num_X + num_zeta1 + num_zeta2 + num_v + num_z);
    
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

    % ---- Construct zprev cost Constraint: -z ≤ -z_prev ----
    A_z = zeros(num_z,num_U + num_X + num_zeta1 + num_zeta2 + 1 + num_z);
    A_z(end) = -1;
    b_z = -z_prev;
    % ---- Construct terminal cost Constraint: -v + a_i*z ≤ -b_i ----
    
    % Initialize terminal cost constraint matrix 
    A_term = zeros(num_i,num_U + num_X + num_zeta1 + num_zeta2 + 1 + num_z);
    
    % Assign -1 vector
    A_term(:,end-1:end-1) = -1*ones(num_i,1);
    % Assign a_i vector 
    A_term(:,end) = a_i; 

    % Right-hand side for terminal cost constraint
    b_term = -b_i;

    % ---- Combine Constraints ----
    E_in = [A_temp; A_peak; A_z; A_term];
    e_in = [b_temp; b_peak; b_z; b_term];

end

