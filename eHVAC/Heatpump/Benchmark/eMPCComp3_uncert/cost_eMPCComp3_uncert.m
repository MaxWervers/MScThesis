function bill = cost_eMPCComp3_uncert(Nh_add, start, month, num_days, meta_data)

    % addpath('MetaData')
    % load('MetaData\meta_data.mat','meta_data')
        
    [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump();
    
    
    % Define system parameters
    H_split = 4;                % Split hour into 4 (15 min intervals)
    
    [N, N_vec_known, N_vec_ext, Steps, Ta, c, Tmin_air, Tmax_air] = sim_parameters_Nvec3(Nh_add, num_days, start, H_split, meta_data);
    cp = 5*10e-3;               % Peak demand charge eur/wh 
    
    % Disturbance vector
    d = [Ta ...
        ; zeros(1,Steps+H_split)];
   
    % Constraints
    umin = 0;
    umax = 10000;
    Tmax_w = 45;
    
    %%
    % Detect number of predictions
    T_pred_matrix = meta_data{start, 11}{1};  % Assumes same size for all
    num_preds = size(T_pred_matrix, 1);
    bills_array = zeros(num_preds, 1);

    % Loop over each prediction scenario
    for pred_idx = 1:num_preds
        u_prev_peak = 0;
        X_sim = zeros(n, Steps);
        U_sim = zeros(mc, Steps - 1);
        X_sim(:,1) = [20; 20; 30];

        for t = 1:Steps-1

            size_opt = N_vec_ext(t)*(mc+n+2) + 1; % lenght of optimization vector
            E_eq = genE(Ad,Bd,N_vec_ext(t),size_opt);
            
            lb = -10e6*ones(size_opt,1);
            lb(1:N_vec_ext(t)*mc) = umin;
            lb(N_vec_ext(t)*(mc+n)+1:end) = 0; % slack1,2 > 0, z > 0  
            
            ub = 10e6*ones(size_opt,1);
            ub(1:N_vec_ext(t)*mc) = umax;
            
            model.lb = lb;           % Lower bounds
            model.ub = ub;           % Upper bounds
                    

            start_pred_c = start + floor((t + N_vec_known(t) - 1)/H_split);
            c_pred = meta_data{start_pred_c, 10}{1};
            c_pred = repelem(c_pred, H_split);
            c_real_pred = [c(t:t + min(N_vec_known(t), N_vec_ext(t))-1), c_pred(1:max(N_vec_ext(t) - N_vec_known(t), 0))];

            start_pred_T = start + floor(t/H_split);
            T_pred_full = meta_data{start_pred_T, 11}{1};  % 30xN matrix
            T_pred = repelem(T_pred_full(pred_idx, :), H_split);
            d_pred = [T_pred(1:N_vec_ext(t)); d(2:md, t:t+N_vec_ext(t)-1)];

            w = 2;
            a_tel = ((t+N_vec_ext(t)-1)/Steps)^w;

            f = [c_real_pred' * 10e-6 * (1/H_split); zeros(N_vec_ext(t)*n, 1); 10e8 * ones(2*N_vec_ext(t), 1); cp*a_tel];
            e_eq = gen_e(Ad, Bd_d, d_pred, N_vec_ext(t), X_sim(:,t));
            [E_in, e_in] = construct_soft_constraints_HP(N_vec_ext(t), mc, n, Tmin_air(t+1:t+N_vec_ext(t)), Tmax_air(t+1:t+N_vec_ext(t)), Tmax_w, u_prev_peak);

            
            model.obj = f;
            model.A = [sparse(E_in); sparse(E_eq)];
            model.rhs = [e_in; e_eq];
            model.sense = [repmat('<', size(e_in)); repmat('=', size(e_eq))];
            model.modelsense = 'min';

            params.outputflag = 0;
            result = gurobi(model, params);
            sol = result.x;

            U_sim(:,t) = sol(1);
            u_prev_peak = max(u_prev_peak,sol(1));
            X_sim(:,t+1) = Ad * X_sim(:,t) + Bd * U_sim(:,t) + Bd_d * d(:,t);
        end

        % Compute bill for this prediction
        bills_array(pred_idx) = dot(c(1:length(U_sim)), U_sim) * (1/H_split) * 10e-6 + max(U_sim) * cp;
    end

    bill = bills_array;

    vars = setdiff(who, 'meta_data');
    save(['Heatpump\Benchmark\eMPCComp3_uncert\cost_eMPCComp3_uncert',num2str(Nh_add),'h_add_', month, '_FULL'], vars{:})
    
 end