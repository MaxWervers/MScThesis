function bill = optimal_cost_eMPC5(num_days, MONTH, start, meta_data)

    % addpath('MetaData')
    % load('MetaData\meta_data.mat','meta_data')

    % Dynamic model
    [Ad, Bd, Bd_d, n, mc, md] = dyn_heatpump();

    % Time resolution setup
    H_split = 4;                             % 15-min intervals
    Nh = num_days * 24;                      % Hours

    % Get simulation parameters
    [N, Steps, Ta, c, Tmin_air, Tmax_air] = sim_parameters(Nh, num_days, start, H_split, meta_data);
    cp = 5*10e-3;                            % Peak demand charge eur/Wh

    % Disturbance vector
    d = [Ta; zeros(1, Steps + N)];

    % Initial state
    X_sim = zeros(n, N+1);
    U_sim = zeros(mc, N);
    X_sim(:,1) = [20; 20; 30];

    % Constraints
    umin = 0;
    umax = 20000;
    Tmax_w = 45;

    % Optimization setup
    size_opt = N * (mc + n + 2) + 1;
    E_eq = genE(Ad, Bd, N, size_opt);

    lb = -10e6 * ones(size_opt,1);
    ub =  10e6 * ones(size_opt,1);
    lb(1:N*mc) = umin;
    ub(1:N*mc) = umax;
    lb(N*(mc+n)+1:end) = 0;

    model.lb = lb;
    model.ub = ub;

    % Cost vector
    f = [c(1:N)' * 10e-6 / H_split; zeros(N*n,1); 10e8 * ones(2*N,1); cp];

    % Constraints
    e_eq = gen_e(Ad, Bd_d, d(1:md, 1:N), N, X_sim(:,1));
    [E_in, e_in] = construct_soft_constraints_HP(N, mc, n, Tmin_air(2:N+1), Tmax_air(2:N+1), Tmax_w, 0);

    % Set up Gurobi model
    model.obj = f;
    model.A = [sparse(E_in); sparse(E_eq)];
    model.rhs = [e_in; e_eq];
    model.sense = [repmat('<', size(e_in)); repmat('=', size(e_eq))];
    model.modelsense = 'min';

    % Gurobi parameters
    params.outputflag = 0;

    % Solve optimization
    result = gurobi(model, params);
    sol = result.x;

    % Extract result
    U_sim(1:N) = sol(1:N);
    X_sim(:, 2:N+1) = reshape(sol(N+1:N+3*N), n, []);
    z = sol(end);

    % Compute bill
    bill = dot(c(1:N), U_sim) * (1/H_split) * 10e-6 + max(U_sim) * cp;

    vars = setdiff(who, 'meta_data');
    save(['Heatpump\Benchmark\OptimalCost\opti_cost_HPeMPC5_', MONTH, '_FULL'], vars{:})

end

