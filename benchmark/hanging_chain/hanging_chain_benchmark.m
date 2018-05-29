clear all;

addpath('./model/');

N = 40;
Ts = 0.2;
num_sim_iters = 50;

mass_vs_time = [];

for num_masses=3:3


    num_free_masses = num_masses - 2;
    p_end_ref = [1; 0; 0];
    v_ref = zeros(3*num_free_masses, 1);

    [p_ref, ~, ~] = fsolve(@(p) hanging_chain_ode([p; v_ref; p_end_ref], zeros(3, 1), num_free_masses), ...
                    linspace(0, 1, 3*num_free_masses).', ...   
                    optimoptions('fsolve','Algorithm','Levenberg-Marquardt'));

                
    x_ref = [p_ref; v_ref; p_end_ref];

    x0 = x_ref;
    x0(end-2:end) = [1; -1; 0];

    Q = blkdiag(25*eye(length(p_ref)), 1*eye(length(v_ref)), 25*eye(3));
    R = 0.01*eye(3);

    W = blkdiag(Q, R);
    WN = Q;

    umax = 1;

    solvers = {'acados', 'acado'};

    integrator = @(x, u) ode45(@(t, y) hanging_chain_ode(y, u, num_free_masses), [0 Ts], x);

    for num_solver=1:numel(solvers)
        sol = solvers{num_solver};

        [~,~,~] = mkdir(['_', sol]);
        copyfile([sol, '_run.m'], ['_', sol]);
        eval(['cd _', sol]);
        eval(['[X.', sol, ', U.', sol, ', timing.', sol, ', status.', sol, ...
            ', num_iters.', sol, ']', ...
            '= ', sol, '_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x0, num_sim_iters, integrator);']);
        cd ..

        mass_vs_time = [mass_vs_time; mean(timing.(sol))];

    end


end

mass_vs_time = reshape(mass_vs_time, numel(solvers), []);

analysis;
