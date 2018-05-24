
clear all;

addpath('model');

N = 20;

nx = 4;
nu = 1;

Ts = 0.05;
Q = diag([1, 1, 1e-10, 1e-10]);
W = blkdiag(Q, 1e-2);
WN = 1000*Q;

Fmax = 8;
x0 = [0; pi; 0; 0];
num_sim_iters = 100;

solvers = {'acados', 'falcopt', 'grampc', 'ipopt', 'acado'};

for num_solver=1:numel(solvers)
    sol = solvers{num_solver};

    [~,~,~] = mkdir(['_', sol]);
    copyfile([sol, '_run.m'], ['_', sol]);
    eval(['cd _', sol]);
    eval(['[X.', sol, ', U.', sol, ', timing.', sol, ', status.', sol, ']', ...
          '= ', sol, '_run(N, Ts, W, WN, Fmax, x0, num_sim_iters);']);
    cd ..
end

analysis;
