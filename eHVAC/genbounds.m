function [lb, ub] = genbounds(N, Bd, Tmin, Tmax, umin, umax)
    % Get dimensions
    [n, mc] = size(Bd);
    
    % Total decision variable size (N*m + N*n)
    num_u = N * mc;
    num_x = N * n;
    
    % Initialize bounds with large values
    lb = -1e5 * ones(num_u + num_x, 1);
    ub = 1e5 * ones(num_u + num_x, 1);
    
    % Set bounds for control inputs
    for k = 1:N
        row_start = (k-1) * mc + 1;
        row_end = row_start + mc - 1;
        lb(row_start:row_end) = umin;
        ub(row_start:row_end) = umax;
    end
    
    % Set bounds for the second state variable
    for k = 1:N
        row_start = num_u + (k-1) * n + 2; % +2 enforces second state index
        lb(row_start) = Tmin(k);
        ub(row_start) = Tmax(k);
    end
end

