function e = gen_e(Ad, Bd_d, d_vec, N, x0)
    % GEN_E - Generate e vector for the optimization problem
    %
    % Syntax:  e = gen_e(Ad, Bd_d, d_vec, N, x0)
    %
    % Inputs:
    %    Ad - Discrete-time state matrix n x n
    %    Bd_d - Discrete-time disturbance input matrix n x md
    %    d_vec - Disturbance vector md x N
    %    N - Prediction horizon, in steps (scalar)
    %    x0 - Initial state vector n x 1
    %
    % Outputs:
    %    e - e vector for the optimization problem

    % Get dimensions
    [n, md] = size(Bd_d);
    
    % Initialize e vector with Bd_ds
    e = zeros(N * n, 1);
    
    % Fill e vector row by row for k = 0 to N-1
    for k = 1:N
        row_start = (k-1) * n + 1;
        row_end = row_start + n - 1;
        
        % Compute e vector using d(k)
        e(row_start:row_end) = Bd_d*d_vec(:,k);
    end
    
    % Adjust first block of e to include A_d * x0
    e(1:n) = e(1:n) + Ad * x0;


end