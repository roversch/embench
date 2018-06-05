function [states, controls, timing, status, num_iters] = acado_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

rng(0);

clear global;

eval(['DifferentialState p(',num2str(num_free_masses*3),',1);']);  % 3-dimensional position of masses 1, 2, ..., M
DifferentialState p_end(3,1);  % 3-dimensional position of end point (M+1)
eval(['DifferentialState v(',num2str(num_free_masses*3),',1);']);  % 3-dimensional velocity of masses 1, 2, ..., M
Control u(3,1);

% Differential Equation

ode_rhs = hanging_chain_ode([p; p_end; v], u, num_free_masses);

acadoSet('problemname', 'mpc');

ocp = acado.OCP( 0.0, N*Ts, N );

ocp.minimizeLSQ( W, [p; p_end; v; u] );
ocp.minimizeLSQEndTerm( WN, [p; p_end; v] );

ocp.subjectTo( -umax <= u <= umax ); % bounds on controls

ocp.setModel(ode_rhs);

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',        N                   );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);

mpc.exportCode('chain_export');

cd('chain_export');

global ACADO_;
copyfile([ACADO_.pwd '/../../external_packages/qpoases'], './qpoases');
copyfile([ACADO_.pwd '/../../external_packages/qpoases3'], './qpoases3');

make_acado_solver('../acado_MPCstep');
cd('../');

avg_timing = 0;

for irepeat=1:nrepeat

input.x = repmat(x_0.', N+1, 1);
input.u = repmat(zeros(1, 3), N, 1);
input.y = repmat([x_ref.', zeros(1, 3)], N, 1);
input.yN = x_ref.';

states = x_0.';
controls = [];
timing = [];
status = [];
num_iters = [];

for i=1:num_sim_iters
    
    input.x0 = states(end, :);
    
    total_time = 0;
    KKT=1e8;
    KKTiter = 0;
    while KKT > 1e-9 && KKTiter < 1
        output = acado_MPCstep(input);
        input.x = output.x;
        input.u = output.u;
        KKT = output.info.kktValue;
        total_time = total_time + output.info.cpuTime;
        KKTiter = KKTiter + 1;
    end
    
    input.x = output.x;
    input.u = output.u;
    
    if i >= num_sim_iters/2 && i < num_sim_iters/2 + 5
        controls = [controls; u_0.'];
    else
        controls = [controls; output.u(1, :)];
    end
    
    [~, sim_out] = integrator_fun(input.x0.', controls(end, :).');
    
    states = [states; sim_out(end, :) + sigma*randn(1, 3*(2*num_free_masses+1))];
    timing = [timing; total_time];
    status = [status; output.info.status];
    num_iters = [num_iters; output.info.QP_iter];
    
end

avg_timing = avg_timing + timing;

end

timing = avg_timing / nrepeat;

system('rm acado_MPCstep.mexmaci64');

end

