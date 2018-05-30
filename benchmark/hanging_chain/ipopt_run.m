function [states, controls, timing, status, num_iters] = ipopt_run(num_free_masses, N, Ts, W, WN, umax, x_ref, u_0, num_sim_iters, integrator_fun)

clear GLOBAL

import casadi.*

addpath('~/local/matlab');
addpath('~/local/lib');

nx = (2*num_free_masses+1)*3;
nu = 3;

x = SX.sym('x', nx, 1);
u = SX.sym('u', nu, 1);

ode = Function('ode', {x, u}, {hanging_chain_ode(x, u, num_free_masses)});
F = simpleRK(ode, 1, 4);

xu_ref = [x_ref; zeros(nu, 1)];

f = 0.0;
x = MX.sym('x_0', nx);
w = x;
w0 = x_ref;
g = [];
for i=0:N-1
    
    x = w(end-nx+1:end);
    u = MX.sym(['u_', num2str(i)], nu);
    x1 = MX.sym(['x_', num2str(i+1)], nx);

    sim_out = F('x0', x, 'p', u, 'h', Ts);
    
    w = [w; u; x1];
    w0 = [w0; zeros(nu, 1); x_ref];
    
    g = [g; x1 - sim_out.xf];
    
    xu = [x; u];
    f = f + (xu - xu_ref).' * W * (xu - xu_ref);
    
end

x = w(end-nx+1:end);
f = f + (x-x_ref).' * WN * (x-x_ref);

ipopt_solver = nlpsol('ipopt_solver', 'ipopt', struct('x', w, 'f', f, 'g', g), ...
                        struct('expand', 1, 'print_time', 0, ...
                        'ipopt', struct('print_level', 1, 'linear_solver', 'mumps')));

bound = [umax; +inf*ones(nx, 1)];

states = x_ref.';
controls = [];
timing = [];
status = [];
num_iters = [];

for i=1:num_sim_iters

    x0 = states(end, :).';
    
    lbx = [x0; repmat(-bound, N, 1)];
    ubx = [x0; repmat(+bound, N, 1)];

    output = ipopt_solver('x0', w0, 'lbx', lbx, 'ubx', ubx, 'lbg', 0, 'ubg', 0);

    w_opt = full(output.x);
    w0 = w_opt;

    if (i <= 5)
        controls = [controls; u_0.'];
    else
        controls = [controls; w_opt(nx+1: nx+nu).'];
    end

    [~, sim_out] = integrator_fun(x0, controls(end, :).');

    states = [states; sim_out(end, :)];

    timing = [timing; ipopt_solver.stats.t_proc_ipopt_solver];

    if (strcmp(ipopt_solver.stats.return_status, 'Solve_Succeeded'))
        status = [status; 0];
    else
        status = [status; -1];
    end

    num_iters = [num_iters; ipopt_solver.stats.iter_count];
    
end

end
    
    