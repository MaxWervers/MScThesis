function cost_to_go_v3(num_days,MONTH,start)
    % Function to simulate heating system costs over varying days and z_fix values
    % Inputs:
    %   num_days - Number of days to simulate
    %   start    - Start index into meta_data

    addpath('MetaData')
    load('MetaData\meta_data.mat','meta_data')

    %% Dynamic model
    [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump();

    %% Parameters
    H_split = 4;                        % 15 min intervals
    end_index = start + num_days*24;         % Keep fixed 30-day end
    z_fix_values = [10:10:140,150:25:475 ,500:50:1000];
    num_days_values = num_days:-1:1;   % Days decreasing from input
    cp = 5*10e-3;  % Peak demand charge eur/wh 

    lag = 5; %3 = 48h nhadd

    % Storage
    results = cell(length(num_days_values),1);  % cell array to store results for each num_days_j
    tol = 1e-6;  % tolerance for detecting no change in bill

    % Loop over each num_days_values first
    for j = 1:length(num_days_values)
        num_days_j = num_days_values(j);  % use a specific num_days value
        current_z_fix_plot = [];    % will store the z_fix values for this num_days_j
        current_u_max = [];  % will store the corresponding sol(end) values
        current_bill = [];     % will store the corresponding bill values
        previous_bill = nan;  % initialize tracking for sol(end) for this num_days
        
        % Loop over the z_fix_values for the given num_days_j
        for i = 1:length(z_fix_values)
            z_fix = z_fix_values(i);
            start_j = end_index - (num_days_j * 24);
            start_j_cpred = end_index - ((num_days_j + lag)* 24);
            Nh = end_index - start_j;
        
            [N, Steps, Ta, c, Tmin_air, Tmax_air, Ta_pred] = ...
                sim_parameters_eMPC8(Nh, num_days_j, start_j, H_split, meta_data, 8);

            % Price prediction
            c_pred = meta_data{start_j_cpred, 10};
            c_pred = repelem(c_pred{1}(1:24), H_split);
            % c_pred = repmat(c_pred, 1, num_days_j-1); %v2
            % c_real_pred = [c(1:24*H_split), c_pred];  %v2
            c_real_pred = repmat(c_pred, 1, num_days_j);%v21

            % Disturbance vector
            d = [Ta_pred(1:N); zeros(1,N)];

            % Initialize simulation
            X_sim = zeros(n, N+1);
            U_sim = zeros(mc, N);
            X_sim(:,1) = [20; 20; 30];

            % Constraints
            umin = 0;
            umax = 20000;
            Tmax_w = 45;
        
            size_opt = N*(mc+n+2) + 1;
            E_eq = genE(Ad, Bd, N, size_opt);
        
            lb = -10e6 * ones(size_opt,1);
            lb(1:N*mc) = umin;
            lb(N*(mc+n)+1:end) = 0;
            lb(end) = z_fix;
        
            ub = 10e6 * ones(size_opt,1);
            ub(1:N*mc) = umax;
            ub(end) = umax;
        
            model.lb = lb;
            model.ub = ub;
        
            % Optimization (single step)
            for t = 1:1
                f = [c_real_pred'*10e-6*(1/H_split); zeros(N*n,1); 10e8*ones(2*N,1); cp];
                e_eq = gen_e(Ad, Bd_d, d(1:md, t:t+N-1), N, X_sim(:,t));
                [E_in, e_in] = construct_soft_constraints_HP(N, mc, n, ...
                    Tmin_air(t+1:t+N), Tmax_air(t+1:t+N), Tmax_w, 0);

                model.obj = f;
                model.A = [sparse(E_in); sparse(E_eq)];
                model.rhs = [e_in; e_eq];
                model.sense = [repmat('<', size(e_in)); repmat('=', size(e_eq))];
                model.modelsense = 'min';

                params.outputflag = 0;
                result = gurobi(model, params);

                sol = result.x;
                U_sim(1:N) = sol(1:N);
                X_sim(:,2:N+1) = reshape(sol(N+1:N+3*N), n, []);
            end

            z = sol(end);
            u_max = max(U_sim);
            cost = result.objval;
            bill_zmax = dot(c_real_pred', U_sim) * (1/H_split) * 10e-6 + z * cp;
            bill = dot(c_real_pred', U_sim) * (1/H_split) * 10e-6 + u_max * cp;
            
            % Check if constraints can be met otherwise go to next z value
            % if abs(cost - bill_zmax) > 1
            %     % fprintf('Breaking inner loop for z = %.1f, num_day = %.1f, cost = %.4f, bill = %.4f\n', z_fix, num_days_j, cost, bill);
            %     continue;  % go to next z_fix value if constraints are not met
            % end
        
            % If bill is zero, recalculate using z_fix
            % if bill == 0
            %     bill = z_fix * cp;
            %     u_max = z_fix;
            % end
            

            % Check if the last control variable remains unchanged relative to previous z_fix
            % if ~isnan(previous_bill)
            %     if abs(bill - previous_bill) < tol
            %         fprintf('bill did not change for new z_fix = %.1f (num_day = %.1f). Breaking further z_fix iterations.\n', z_fix, num_days_j);
            %         break;  % exit the z_fix loop for this num_days_j value
            %     end
            % end
            previous_bill = bill; % update tracker for next iteration

            current_z_fix_plot = [current_z_fix_plot, z_fix];
            current_u_max = [current_u_max, u_max];
            current_bill = [current_bill, cost];

        end
        
        % If no valid z_fix values were stored, fill with defaults using all z_fix_values.
        % if isempty(current_z_fix)
        %     current_z_fix = z_fix_values; 
        %     current_u_max = z_fix_values; % Here we assign u_max = z_fix for all.
        %     current_bill = current_z_fix * cp;  % bill = z_fix * cp element-wise
        % end

        % Store the results for the current num_days_j in the cell array:
        results{j} = struct('num_days', num_days_j, ...
                        'z_fix', current_z_fix_plot, ...
                        'u_max', current_u_max, ...
                        'bill', current_bill); 
    end

% Save results with MONTH in filename
vars = setdiff(who, 'meta_data');
save(['Heatpump\Benchmark\TermCost_v3\cost_to_go_v3_96h_', MONTH, '_FULL'], vars{:})
save(['Heatpump\Benchmark\TermCost_v3\cost_to_go_v3_96h_', MONTH], "results")
end