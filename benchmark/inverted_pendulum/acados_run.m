function [states, controls, timing, status, num_iters] = acados_run(N,...
    Ts, W, WN, Fmax, x0, num_sim_iters)

import casadi.*

model = pendulum_ode_model;
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

% cost
cost_expr = [x;u]'* W * [x;u];
ocp_model.set('cost_expr_ext_cost', cost_expr);
ocp_model.set('cost_expr_ext_cost_e', x' * WN * x );

% dynamics
ocp_model.set('dyn_type', 'explicit');
ocp_model.set('dyn_expr_f', model.expr_f_expl);

% constraints
ocp_model.set('constr_type', 'bgh');
ocp_model.set('constr_Jbu', eye(nu));
ocp_model.set('constr_lbu', -Fmax); % lower bound on u
ocp_model.set('constr_ubu', Fmax);  % upper bound on u

ocp_model.set('constr_x0', x0);

%% acados ocp set optsbuffer
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', N);
ocp_opts.set('nlp_solver', 'sqp_rti');
ocp_opts.set('sim_method', sim_method);
qp_solver = 'partial_condensing_hpipm';
ocp_opts.set('qp_solver', qp_solver);
qp_solver_cond_N = round(N/4); % for partial condensing
ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);

ocp_opts.set('nlp_solver_exact_hessian', 'false');

%% create ocp solver
ocp = acados_ocp(ocp_model, ocp_opts);
for j=1:1

    states = x0.';
    controls = [];
    timing = [];
    status = [];
    num_iters = [];

%     x_traj_init = zeros(nx, N+1);
%     ocp.set('init_x', x_traj_init);

    for i=1:num_sim_iters

        current_x = states(end, :).'
        
        ocp.set('constr_x0', current_x);

        ocp.solve();
        
        xnext = ocp.get('x', 1);
        unext = ocp.get('u', 1);
        states = [states; xnext'];
        controls = [controls; unext];

        timing = [timing; ocp.get('time_tot') ];
        if ocp.get('status')
            keyboard
        end
        status = [status; ocp.get('status')];
        
        stat = ocp.get('stat');
        qp_iter = sum(stat(:,3));
        num_iters = [num_iters; qp_iter];
            
    end

end
        
end

