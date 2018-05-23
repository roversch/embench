function [states, controls, timing, status] = falcopt_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

clear GLOBAL

import casadi.*

dynamics = @(x,u) pendulum_ode_discrete(x, u, Ts);

nx = 4;
nu = 1;

objective = struct('Q', W(1:nx, 1:nx),'R', W(end, end), 'P', WN);

variable_stepSize.active = true;

info = falcopt.generateCode(dynamics, N, nx, nu, objective, ...
    'variable_stepSize', variable_stepSize, ...
    'gradients', 'casadi', ...
    'lb', -Fmax, 'ub', +Fmax, ...
    'contractive', false, 'terminal', false, ...
    'debug', 3, 'merit_function', 0,...
    'eps', 1e-3, 'precision', 'double', ...
    'name', 'pendulum_falcopt', 'gendir', 'pendulum_export');


states = x0.';
controls = [];
timing = [];
status = [];

u_guess = zeros(N, nu);

for i=1:num_sim_iters
    
    [~, flag, info] = pendulum_falcopt(states(end, :).', u_guess);
    
    states = [states; info.x(:, 1).'];
    controls = [controls; info.u(:, 1).'];
    timing = [timing; info.time];
    if flag > 0
        status = [status; 0];
    else
        status = [status; -1];
    end
    
    u_guess = info.u.';

end

end

