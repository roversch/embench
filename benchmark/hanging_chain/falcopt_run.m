function [states, controls, timing, status, num_iters] = falcopt_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

rng(0);

clear GLOBAL

addpath('~/FalcOpt');
addpath('~/local/matlab');

import casadi.*

dynamics = @(x,u) hanging_chain_ode_discrete(x, u, Ts, num_free_masses);

nx = (num_free_masses * 2 + 1) * 3;
nu = 3;

objective = struct('Q', W(1:nx, 1:nx),'R', W(end-nu+1:end, end-nu+1:end), 'P', WN);
objective.trackReference = true;

variable_stepSize.active = false;
% variable_stepSize.alpha_max = 1;
%     variable_stepSize.steady_state_state = x_ref;
%     variable_stepSize.steady_state_input = zeros(3, 1);

info = falcopt.generateCode(dynamics, N, nx, nu, objective, ...
    'variable_stepSize', variable_stepSize, ...
    'gradients', 'casadi', ...
    'lb', -umax, 'ub', +umax, ...
    'contractive', false, 'terminal', false, ...
    'debug', 3, 'merit_function', 0, 'maxIt', 100, ...
    'eps', 1e-1, 'precision', 'double', ...
    'name', 'hanging_chain_falcopt', 'gendir', 'hanging_chain_export');

all_xref = repmat(x_ref.', N, 1).';
all_uref = repmat(zeros(1,3), N, 1).';

avg_timing = 0;

for irepeat=1:nrepeat

states = x_0.';
controls = [];
timing = [];
status = [];
num_iters = [];

u_guess = zeros(N, nu);

for i=1:num_sim_iters

    [~, flag, info] = hanging_chain_falcopt(states(end, :).', all_xref(:), all_uref(:), u_guess(:));
    
    if i >= num_sim_iters/2 && i < num_sim_iters/2 + 5
        controls = [controls; u_0.'];
    else
        controls = [controls; info.u(:, 1).'];
    end
    [~, sim_out] = integrator_fun(states(end, :).', controls(end, :).');

    states = [states; sim_out(end, :) + sigma*randn(1, 3*(2*num_free_masses+1))];
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

avg_timing = avg_timing + timing;

end

timing = avg_timing / nrepeat;

system('rm hanging_chain_falcopt.mexmaci64');
    
end
    
    