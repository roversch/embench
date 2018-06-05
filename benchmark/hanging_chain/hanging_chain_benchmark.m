clear all;

addpath('./model/');

sigma = 0;
nrepeat = 10;

N = 40;
Ts = 0.2;
num_sim_iters = round(60/Ts);

mass_vs_time = [];

u_0 = [-1; 1; 1];

hp = zeros(20, 20);

for num_masses=5:5

    num_free_masses = num_masses - 2;
    p_end_ref = [7.5; 0; 0];
    v_ref = zeros(3*num_free_masses, 1);

    [p_ref, ~, ~] = fsolve(@(p) hanging_chain_ode([p; p_end_ref; v_ref], zeros(3, 1), num_free_masses), ...
                    linspace(0, 1, 3*num_free_masses).', ...   
                    optimoptions('fsolve','Algorithm','Levenberg-Marquardt'));

                
    x_equilibrium = [p_ref; p_end_ref; v_ref];
    x_ref = [0*p_ref; p_end_ref; 0*v_ref];

    x_0 = 0*x_ref;
    for kk = 1:num_free_masses+1
        x_0((kk-1)* 3 + 1) = 7.5 * kk / (num_free_masses+1);
    end
    
    Q = blkdiag(1e-10*eye(length(p_ref)), 2.5*eye(3), 25*eye(length(v_ref)));
    R = 0.1*eye(3);

    W = blkdiag(Q, R);
    WN = blkdiag(1e-10*eye(length(p_ref)), 10*eye(3), 1e-10*eye(length(v_ref)));

    umax = ones(3, 1);

    solvers = {'ipopt', 'falcopt', 'viatoc', 'acado', 'grampc', 'acados'};

    integrator = @(x, u) ode45(@(t, y) hanging_chain_ode(y, u, num_free_masses), [0 Ts], x);

    for num_solver=1:numel(solvers)
        sol = solvers{num_solver};

        [~,~,~] = mkdir(['_', sol]);
        copyfile([sol, '_run.m'], ['_', sol]);
        eval(['cd _', sol]);
        eval(['[X.', sol, ', U.', sol, ', timing.', sol, ', status.', sol, ...
            ', num_iters.', sol, ']', ...
            '= ', sol, '_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator, sigma, nrepeat);']);
        cd ..

        mass_vs_time = [mass_vs_time; mean(timing.(sol))];

    end

end

mass_vs_time = reshape(mass_vs_time, numel(solvers), []);

analysis;
