function [states, controls, timing, status, num_iters]  = acados_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

rng(0);

addpath('~/local/matlab');
addpath('~/local/lib');

import casadi.*
import acados.*

p = SX.sym('p', 3*num_free_masses);
p_end = SX.sym('p_end', 3);
v = SX.sym('v', 3*num_free_masses);
u = SX.sym('u', 3);

ode_rhs = hanging_chain_ode([p; p_end; v], u, num_free_masses);
ode_fun = Function('ode_fun', {[p; p_end; v], u}, {ode_rhs});

nx = (num_free_masses * 2 + 1) * 3;
nu = 3;

nlp = ocp_nlp(N, nx, nu);
nlp.set_dynamics(ode_fun, struct('integrator', 'rk4', 'step', Ts));

nlp.set_field('lbu', -umax);
nlp.set_field('ubu', +umax);

nlp.set_field('lbx', 0, x_0);
nlp.set_field('ubx', 0, x_0);

nlp.set_stage_cost(eye(nx+nu), [x_ref; zeros(nu, 1)], W);
nlp.set_terminal_cost(eye(nx), x_ref, WN);

nlp.initialize_solver('sqp', struct('qp_solver', 'hpipm', 'max_iter', 1, ...
                      'hpipm', struct('N2', 5)));

avg_timing = 0;
                  
for irepeat=1:nrepeat
                  
states = x_0.';
controls = [];
timing = [];
status = [];
num_iters = [];

for i=1:num_sim_iters

    current_x = states(end, :).';

    nlp.set_field('lbx', 0, current_x);
    nlp.set_field('ubx', 0, current_x);
    if (i == 1)
        output = nlp.solve(x_0, zeros(3, 1));
    else
        output = nlp.solve();
    end

    all_x = output.states();
    all_u = output.controls();

    if i >= num_sim_iters/2 && i < num_sim_iters/2 + 5
        controls = [controls; u_0.'];
    else
        controls = [controls; all_u{1}.'];
    end

    [~, sim_out] = integrator_fun(current_x, controls(end, :).');
    states = [states; sim_out(end, :) + sigma*randn(1, 3*(2*num_free_masses+1))];
    
    info = output.info();

    timing = [timing; info.total_time];
    status = [status; info.status];

    num_iters = [num_iters; info.num_qp_iter];

end

avg_timing = avg_timing + timing;

end

timing = avg_timing / nrepeat;

end

