function E = genE(Ad, Bd, N, size_opt)
    % GEN_E - Generate E matrix for the optimization problem
    %
    % Syntax:  E = genE(Ad, Bd, N)
    %
    % Inputs:
    %    Ad - Discrete-time state matrix
    %    Bd - Discrete-time control input matrix
    %    N - Prediction horizon, in steps
    %
    % Outputs:
    %    E - E matrix for the optimization problem
    %
    
    % Get dimensions
    [n, mc] = size(Bd);
    
    % Initialize E matrix of size (N*n) x (N*m + N*n)
    E = zeros(N * n, size_opt);
    
    % Fill E matrix row by row for k = 0 to N-1
    for k = 1:N
        row_start = (k-1) * n + 1;
        row_end = row_start + n - 1;
        
        % Control input block (-Bd)
        col_u_start = (k-1) * mc + 1;
        col_u_end = col_u_start + mc - 1;
        E(row_start:row_end, col_u_start:col_u_end) = -Bd;
        
        % State variable block
        col_x_start = N*mc + (k-1) * n + 1;
        col_x_end = col_x_start + n - 1;
        
        % Current state (-Ad)
        if k > 1
            prev_x_start = N*mc + (k-2) * n + 1;
            prev_x_end = prev_x_start + n - 1;
            E(row_start:row_end, prev_x_start:prev_x_end) = -Ad;
        end
        
        % Next state (I_n)
        E(row_start:row_end, col_x_start:col_x_end) = eye(n);

        
    end
end

