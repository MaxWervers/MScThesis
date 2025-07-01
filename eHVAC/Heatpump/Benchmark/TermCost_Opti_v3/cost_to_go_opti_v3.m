function cost_to_go_opti_v3(num_days,MONTH,start)
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

    % Storage
    results = cell(length(num_days_values),1);  % cell array to store results for each num_days_j
   
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
            Nh = end_index - start_j;
        
            [N, Steps, Ta, c, Tmin_air, Tmax_air] = sim_parameters(Nh, num_days_j, start_j, H_split, meta_data);
            cp = 5*10e-3;
        
            d = [Ta(1:N); zeros(1,N)];
        
            X_sim = zeros(n, N+1);
            U_sim = zeros(mc, N);
            X_sim(:,1) = [20; 20; 30];
        
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
                f = [c(t:t+N-1)'*10e-6*(1/H_split); zeros(N*n,1); 10e8*ones(2*N,1); cp];
                e_eq = gen_e(Ad, Bd_d, d(1:md,t:t+N-1), N, X_sim(:,t));
                [E_in, e_in] = construct_soft_constraints_HP(N, mc, n, Tmin_air(t+1:t+N), Tmax_air(t+1:t+N), Tmax_w, 0);
        
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
            bill_zfix = dot(c(1:length(U_sim)), U_sim) * (1/H_split) * 10e-6 + z * cp;
            bill = dot(c(1:length(U_sim)), U_sim) * (1/H_split) * 10e-6 + u_max * cp;
            
           

            previous_bill = bill; % update tracker for next iteration

            current_z_fix_plot = [current_z_fix_plot, z_fix];
            current_u_max = [current_u_max, u_max];
            current_bill = [current_bill, bill_zfix];

        end
        
     

        % Store the results for the current num_days_j in the cell array:
        results{j} = struct('num_days', num_days_j, ...
                        'z_fix', current_z_fix_plot, ...
                        'u_max', current_u_max, ...
                        'bill', current_bill); 
    end

% Save results with MONTH in filename
vars = setdiff(who, 'meta_data');
save(['Heatpump\Benchmark\TermCost_Opti_v3\cost_to_go_opti_v3_', MONTH, '_FULL'], vars{:})
save(['Heatpump\Benchmark\TermCost_Opti_v3\cost_to_go_opti_v3_', MONTH], "results")
end