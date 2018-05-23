function [states, controls, timing, status] = ipopt_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

clear GLOBAL

import casadi.*

nx = 4;
nu = 1;

x = SX.sym('x', nx, 1);
u = SX.sym('u', nu, 1);

ode = Function('ode', {x, u}, {pendulum_ode(x, u)});
F = simpleRK(ode, 1, 4);

f = 0.0;
x = MX.sym('x_0', nx);
w = x;
w0 = x0;
g = [];
for i=0:N-1
    
    x = w(end-nx+1:end);
    u = MX.sym(['u_', num2str(i)], nu);
    x1 = MX.sym(['x_', num2str(i+1)], nx);

    sim_out = F('x0', x, 'p', u, 'h', Ts);
    
    w = [w; u; x1];
    w0 = [w0; 0; x0];
    
    g = [g; x1 - sim_out.xf];
    
    xu = [x; u];
    f = f + xu.' * W * xu;
    
end

x = w(end-nx+1:end);
f = f + x.' * WN * x;

ipopt_solver = nlpsol('ipopt_solver', 'ipopt', struct('x', w, 'f', f, 'g', g), ...
                      struct('expand', 1, 'print_time', 0, ...
                      'ipopt', struct('print_level', 1)));

bound = [Fmax; +inf*ones(nx, 1)];

states = x0.';
controls = [];
timing = [];
status = [];

for i=1:num_sim_iters

    x0 = states(end, :).';
    
    lbx = [x0; repmat(-bound, N, 1)];
    ubx = [x0; repmat(+bound, N, 1)];

    output = ipopt_solver('x0', w0, 'lbx', lbx, 'ubx', ubx, 'lbg', 0, 'ubg', 0);

    w_opt = full(output.x);
    w0 = w_opt;
                
    states = [states; w_opt(nx+nu+1:2*nx+nu).'];
    controls = [controls; w_opt(nx+1: nx+nu)];

    timing = [timing; ipopt_solver.stats.t_proc_ipopt_solver];

    if (strcmp(ipopt_solver.stats.return_status, 'Solve_Succeeded'))
        status = [status; 0];
    else
        status = [status; -1];
    end
    
end

end

