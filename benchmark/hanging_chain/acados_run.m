function [states, controls, timing, status, num_iters]  = acados_run(...
    num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

rng(0);

import casadi.*

model = hanging_chain_ode_model(num_free_masses);
x = model.sym_x;
u = model.sym_u;

nx = length(x);
nu = length(u);

%% model to create the solver
ocp_model = acados_ocp_model();
model_name = 'pendulum';
sim_method = 'erk';

%% acados ocp model
ocp_model.set('name', model_name);
T = Ts * N;
ocp_model.set('T', T);

% symbolics
ocp_model.set('sym_x', x);
ocp_model.set('sym_u', u);

% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.expr_f_expl);

% constraints
ocp_model.set('constr_type', 'bgh');
ocp_model.set('constr_Jbu', eye(nu));
ocp_model.set('constr_lbu', -umax);
ocp_model.set('constr_ubu', umax);

ocp_model.set('constr_x0', x_0);


% nlp.set_stage_cost(eye(nx+nu), [x_ref; zeros(nu, 1)], W);
% nlp.set_terminal_cost(eye(nx), x_ref, WN);

% cost
cost_expr = [x-x_ref;u]'* W * [x-x_ref;u];
ocp_model.set('cost_expr_ext_cost', cost_expr);
ocp_model.set('cost_expr_ext_cost_e', (x-x_ref)' * WN * (x-x_ref));

%% acados ocp set opts
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', 'sqp_rti');
ocp_opts.set('sim_method', sim_method);
qp_solver = 'partial_condensing_hpipm';
ocp_opts.set('qp_solver', qp_solver);
qp_solver_cond_N = round(N/8); % for partial condensing
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);

ocp_opts.set('nlp_solver_exact_hessian', 'false');

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);

avg_timing = 0;
                  
for irepeat=1:nrepeat
                  
states = x_0.';
controls = [];
timing = [];
status = [];
num_iters = [];

for i=1:num_sim_iters

    current_x = states(end, :).';

    ocp.set('constr_x0', current_x);
    if (i == 1)
        ocp.set('init_u', repmat(u_0, N, 1));
        ocp.set('init_x', repmat(x_0, N+1, 1));
    end
    ocp.solve();

    u_ocp = ocp.get('u',0);

    % jonathan: why this if?!
    if i >= num_sim_iters/2 && i < num_sim_iters/2 + 5
        controls = [controls; u_0.'];
    else
        controls = [controls; u_ocp.'];
    end

    [~, sim_out] = integrator_fun(current_x, controls(end, :).');
    
    % jonathan: randomness should be the same for all solvers!
    states = [states; sim_out(end, :) + sigma*...
                randn(1, 3*(2*num_free_masses+1))];

    timing = [timing; ocp.get('time_tot') ];
    if ocp.get('status')
        keyboard
    end
    status = [status; ocp.get('status')];

    stat = ocp.get('stat');
    qp_iter = sum(stat(:,3));
    num_iters = [num_iters; qp_iter];

end

avg_timing = avg_timing + timing;

end

timing = avg_timing / nrepeat;

end

