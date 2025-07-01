function bill = cost_eMPCComp2_cert(Nh_add, start, month, num_days, meta_data)

    % addpath('MetaData')
    % load('MetaData\meta_data.mat','meta_data')
    
    %% Dynamic model, discretized with ZOH
    
    [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump();
    
    %% Configure Optimization parameters
   
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
    
    

    %%
    u_prev_peak = 0;
    N_Vec = N_vec_ext;    
    
    for t = 1:Steps-1
    % ========== DEFINE OPTIMIZATION PROBLEM AND RUN MPC LOOP ==========
    
        size_opt = N_Vec(t)*(mc+n+2) + 1; % lenght of optimization vector
        E_eq = genE(Ad,Bd,N_Vec(t),size_opt);
        
        lb = -10e6*ones(size_opt,1);
        lb(1:N_Vec(t)*mc) = umin;
        lb(N_Vec(t)*(mc+n)+1:end) = 0; % slack1,2 > 0, z > 0  
        
        ub = 10e6*ones(size_opt,1);
        ub(1:N_Vec(t)*mc) = umax;
        
        model.lb = lb;           % Lower bounds
        model.ub = ub;           % Upper bounds
    
        sigt = min(((Steps-t)/N),1);
        
        f = [c(t:t+N_Vec(t)-1)'*10e-6*(1/H_split) ; zeros(N_Vec(t)*n,1); 10e8*ones(2*N_Vec(t),1); cp/sigt];
        e_eq = gen_e(Ad, Bd_d, d(1:md,t:t+N_Vec(t)-1), N_Vec(t), X_sim(:,t));
        [E_in, e_in] = construct_soft_constraints_HP(N_Vec(t), mc, n, Tmin_air(t+1:t+N_Vec(t)), Tmax_air(t+1:t+N_Vec(t)), Tmax_w, u_prev_peak);
        
        % f = [c(t:t+N_Vec(t)-1)'*10e-6*(1/H_split); zeros(N_Vec(t)*n,1); 10e8*ones(2*N_Vec(t),1); 1; 0];
        % 
        % e_eq = gen_e(Ad, Bd_d, d(1:md, t:t+N_Vec(t)-1), N_Vec(t), X_sim(:,t));
        % [E_in, e_in] = construct_soft_constraints_HPeMPC9( ...
        %     N_Vec(t), mc, n, Tmin_air(t+1:t+N_Vec(t)), Tmax_air(t+1:t+N_Vec(t)), ...
        %     Tmax_w, 0, u_prev_peak, a_i{tau}', b_i{tau}');
    
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
    

    end

    
    bill = dot(c(1:length(U_sim)),U_sim)*(1/H_split)*10e-6 + max(U_sim)*cp;
    vars = setdiff(who, 'meta_data');
    save(['Heatpump\Benchmark\eMPCComp2_cert\cost_eMPCComp2_cert_',num2str(Nh_add),'h_', month, '_FULL'], vars{:}) 
    
 end




