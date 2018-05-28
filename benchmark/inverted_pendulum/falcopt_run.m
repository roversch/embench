function [states, controls, timing, status, num_iters] = falcopt_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

clear GLOBAL

addpath('~/FalcOpt');

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
num_iters = [];

u_guess = zeros(N, nu);

for i=1:num_sim_iters
    
    [~, flag, info] = pendulum_falcopt(states(end, :).', u_guess);
    
    states = [states; info.x(:, 1).'];
    controls = [controls; info.u(:, 1).'];
    timing = [timing; info.time];
    
    if flag <= 0
        switch flag
            case 0
                if ~isfield(info, 'iterations')
                    warnStr = 'Maximum number of iterations reached. May not be a problem (set ''debug'' to 3 in generateMotor() for more detailed diagnostics). Continuing... ';
                elseif ~isfield(info, 'optimval')
                    warnStr = ['Maximum number of iterations reached with ' num2str(info.iterations) '. May not be a problem (set ''debug'' to 3 in generateMotor() for more detailed diagnostics). Continuing... '];
                else
                    warnStr = ['Maximum number of iterations reached with ' num2str(info.iterations) '. Optimality tolerance ' num2str(info.optimval) ', feasibility tolerance ' num2str(info.feasval) '. Continuing... '];
                end
                warning('falcopt:example:maximumIterations', warnStr);
            case -1
                if ~isfield(info, 'optimval')
                    warnStr = 'Line-search failed to progress. Close to solution, may not be a problem (set ''debug'' to 3 in generateMotor() for more detailed diagnostics). Continuing... ';
                else
                    warnStr = ['Line-search failed to progress. Close to solution, may not be a problem: optimality tolerance ' num2str(info.optimval) ', feasibility tolerance ' num2str(info.feasval) '. Continuing... '];
                end
                warning('falcopt:example:closeToSolution', warnStr);
            case -10
                if ~isfield(info, 'optimval')
                    warnStr = 'Line-search failed to progress (set ''debug'' to 3 in generateMotor() for more detailed diagnostics). Continuing... ';
                else
                    warnStr = ['Line-search failed to progress: optimality tolerance ' num2str(info.optimval) ', feasibility tolerance ' num2str(info.feasval) '. Continuing... '];
                end
                warning('falcopt:example:lineSearchFailed', warnStr);
        end
    end
   
    status = [status; flag - 1];
    
    num_iters = [num_iters; info.iterations];

    u_guess = info.u.';

end

end

