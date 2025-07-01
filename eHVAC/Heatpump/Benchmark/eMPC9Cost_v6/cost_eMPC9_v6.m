function bill = cost_eMPC9_v6(Nh_add, start, month, num_days, meta_data)

    % addpath('MetaData')
    % load('MetaData\meta_data.mat','meta_data')
    % tic
    
    %% Dynamic model, discretized with ZOH
    
    [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump();
    
    %% Configure Optimization parameters
    
    % Define system parameters
    % num_days = 31;              % Number of simulation days
    % Nh_add = 24;
    % start = 105200;             % Start index related to start date in meta_data
    H_split = 4;                % Split hour into 4 (15 min intervals)
    
    [N, N_vec_known, N_vec_ext, Steps, Ta, c, Tmin_air, Tmax_air] = sim_parameters_Nvec3(Nh_add, num_days, start, H_split, meta_data);
    cp = 5*10e-3;               % Peak demand charge eur/wh 
    
    % Disturbance vector
    d = [Ta ...
        ; zeros(1,Steps+H_split)];
    
    % Initialize simulation storage
    X_sim = zeros(n, Steps);
    U_sim = zeros(mc, Steps-1);
    
    % Initial state
    X_sim(:,1) = [20; 20; 30];
    
    % Constraints
    umin = 0;
    umax = 10000;
    Tmax_w = 45;
    
    
    
    %% Terminal cost parameters
    % month = '01-2023';
    filename = ['cost_to_go_opti_v3_', month, '.mat'];
    file_content = open(filename);
    results = file_content.results;
    
    num_days_size = length(results);
    
    % Preallocate as cell arrays
    a_i = cell(num_days_size, 1);
    b_i = cell(num_days_size, 1);
    
    % Loop through each day
    for day = 1:num_days_size
        current_data = results{day};
        
        bill_values = current_data.bill;    % Vector of bill values
        z_values = current_data.z_fix;      % Vector of z_fix values
        
        num_points = length(z_values);
        
        a_day = zeros(1, num_points - 1);
        b_day = zeros(1, num_points - 1);
        
        for i = 1:num_points-1
            a_day(i) = (bill_values(i+1) - bill_values(i)) / ...
                       (z_values(i+1) - z_values(i));
            b_day(i) = bill_values(i) - a_day(i) * z_values(i);
        end
        
        a_i{day} = a_day;
        b_i{day} = b_day;
    end
    
    %%
    % last_printed_progress = -1;
    u_prev_peak = 0;
    % z_prev = 0;
    N_Vec = N_vec_ext;    
    
    for t = 1:Steps-1
    % ========== DEFINE OPTIMIZATION PROBLEM AND RUN MPC LOOP ==========
    
        size_opt = N_Vec(t)*(mc+n+2) + 2; % lenght of optimization vector
        E_eq = genE(Ad,Bd,N_Vec(t),size_opt);
        
        lb = -10e6*ones(size_opt,1);
        lb(1:N_Vec(t)*mc) = umin;
        lb(N_Vec(t)*(mc+n)+1:end) = 0; % slack1,2 > 0, z > 0  
        
        ub = 10e6*ones(size_opt,1);
        ub(1:N_Vec(t)*mc) = umax;
        
        model.lb = lb;           % Lower bounds
        model.ub = ub;           % Upper bounds
    
        % disp((t-1+N_Vec(t))/(H_split*24))
        tau = min(1 + (t-1+N_Vec(t))/(H_split*24), num_days);        
    
        f = [c(t:t+N_Vec(t)-1)'*10e-6*(1/H_split); zeros(N_Vec(t)*n,1); 10e8*ones(2*N_Vec(t),1); 1; 0];
    
        e_eq = gen_e(Ad, Bd_d, d(1:md, t:t+N_Vec(t)-1), N_Vec(t), X_sim(:,t));
        [E_in, e_in] = construct_soft_constraints_HPeMPC9( ...
            N_Vec(t), mc, n, Tmin_air(t+1:t+N_Vec(t)), Tmax_air(t+1:t+N_Vec(t)), ...
            Tmax_w, 0, u_prev_peak, a_i{tau}', b_i{tau}');
    
        model.obj = f;
        model.A = [sparse(E_in); sparse(E_eq)];
        model.rhs = [e_in; e_eq];
        model.sense = [repmat('<', size(e_in)); repmat('=', size(e_eq))];
        model.modelsense = 'min';
    
        params.outputflag = 0;
        result = gurobi(model, params);
    
        sol = result.x;
    
        U_sim(:,t) = sol(1);
        u_prev_peak = max(u_prev_peak, U_sim(:,t));
    
        X_sim(:,t+1) = Ad * X_sim(:,t) + Bd * U_sim(:,t) + Bd_d * d(:,t);
    
    
        % Calculate progress percentage
        % progress = floor((t / (Steps - 1)) * 100); 
    
        % Print only when a new percentage is reached
        % if progress > last_printed_progress
        %     fprintf('Progress: %d%%\n', progress);
        %     last_printed_progress = progress; % Update last printed value
        % end
    end
    % toc
    
    bill = dot(c(1:length(U_sim)),U_sim)*(1/H_split)*10e-6 + max(U_sim)*cp;
    vars = setdiff(who, 'meta_data');
    save(['Heatpump\Benchmark\eMPC9Cost_v6\cost_HPeMPC9_v6_',num2str(Nh_add),'h_add_', month, '_FULL'], vars{:}) %CHANGE BACK
    % %%
    % date_vec_hourly = meta_data.DateTime(start : start+num_days*24-1);
    % minute_offsets = minutes((0:H_split-1) * (60/H_split));  
    % date_vec = repelem(date_vec_hourly, H_split) + repmat(minute_offsets', numel(date_vec_hourly), 1);
    % 
    % figure;
    % subplot(2,1,1);
    % plot(date_vec, X_sim(1,:), 'b', 'LineWidth', 1.5); hold on;
    % plot(date_vec, X_sim(2,:), 'r', 'LineWidth', 1.5);
    % plot(date_vec, X_sim(3,:), 'm', 'LineWidth', 1.5);
    % plot(date_vec, Ta(1:Steps),'g', 'LineWidth', 1);
    % plot(date_vec, Tmin_air(1:Steps), 'Color', [0 0 0], 'LineStyle', '--');
    % plot(date_vec, Tmax_air(1:Steps), 'Color', [0 0 0], 'LineStyle', '--');
    % stairs(date_vec, N_Vec(1:Steps)/4);
    % legend('Tair', 'Tfloor', 'Twater', 'Tamb');
    % xlabel('Time Step');
    % ylabel('State');
    % grid on
    % 
    % subplot(2,1,2);
    % yyaxis left
    % stairs(date_vec(1:end-1), U_sim, 'k', 'LineWidth', 1.5);
    % ylabel('Control Input');
    % yyaxis right
    % stairs(date_vec(1:end-1), c(1:Steps-1), 'r', 'LineWidth', 1.5);
    % ylabel('Elec price');
    % xlabel('Time Step');
 end




