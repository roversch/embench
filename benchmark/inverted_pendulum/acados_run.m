function [states, controls, timing, status] = acados_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

addpath('/Users/robin/local/lib/');

import casadi.*
import acados.*

nx = 4;
nu = 1;

x = SX.sym('x', nx);
u = SX.sym('u', nu);
ode_fun = Function('ode_fun', {x, u}, {pendulum_ode(x, u)});

nlp = ocp_nlp(N, nx, nu);
nlp.set_dynamics(ode_fun, struct('integrator', 'rk4', 'step', Ts));

nlp.set_field('lbu', -Fmax);
nlp.set_field('ubu', +Fmax);

nlp.set_field('lbx', 0, x0);
nlp.set_field('ubx', 0, x0);

nlp.set_stage_cost(eye(nx+nu), zeros(nx+nu, 1), W);
nlp.set_terminal_cost(eye(nx), zeros(nx, 1), WN);

nlp.initialize_solver('sqp', struct('qp_solver', 'hpipm'));

states = x0.';
controls = [];
timing = [];
status = [];

for i=1:num_sim_iters

    current_x = states(end, :).';
    
    nlp.set_field('lbx', 0, current_x);
    nlp.set_field('ubx', 0, current_x);
    if (i == 1)
        output = nlp.solve(x0, 0);
    else
        output = nlp.solve();
    end
    
    all_x = output.states();
    all_u = output.controls();
    states = [states; all_x{2}.'];
    controls = [controls; all_u{1}];
    
    info = output.info();
    
    timing = [timing; info.total_time];
    status = [status; info.status];
    
end
    
end

