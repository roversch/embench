function [states, controls, timing, status, num_iters] = acado_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x0, num_sim_iters, integrator_fun)

clear global;

eval(['DifferentialState p(',num2str(num_free_masses*3),',1);']);  % 3-dimensional position of masses 1, 2, ..., M
eval(['DifferentialState v(',num2str(num_free_masses*3),',1);']);  % 3-dimensional velocity of masses 1, 2, ..., M
DifferentialState p_end(3,1);  % 3-dimensional position of end point (M+1)
Control u(3,1);

% Differential Equation

ode_rhs = hanging_chain_ode([p; v; p_end], u, num_free_masses);

acadoSet('problemname', 'mpc');

ocp = acado.OCP( 0.0, N*Ts, N );

ocp.minimizeLSQ( W, [p; v; p_end; u] );
ocp.minimizeLSQEndTerm( WN, [p; v; p_end] );

ocp.subjectTo( -umax <= u <= umax ); % bounds on controls

ocp.setModel(ode_rhs);

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',        N                   );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES3'    	);

mpc.exportCode('chain_export');

cd('chain_export');

global ACADO_;
copyfile([ACADO_.pwd '/../../external_packages/qpoases'], './qpoases');
copyfile([ACADO_.pwd '/../../external_packages/qpoases3'], './qpoases3');

make_acado_solver('../acado_MPCstep');
cd('../');

input.x = repmat(x_ref.', N+1, 1);
input.u = repmat(zeros(1, 3), N, 1);
input.y = repmat([x_ref.', zeros(1, 3)], N, 1);
input.yN = x_ref.';

states = x0.';
controls = [];
timing = [];
status = [];
num_iters = [];

for i=1:num_sim_iters
    
    input.x0 = states(end, :);
    
    output = acado_MPCstep(input);
    
    input.x = output.x;
    input.u = output.u;
        
    controls = [controls; output.u(1, :)];
    [~, sim_out] = integrator_fun(input.x0.', controls(end, :).');
    
    states = [states; sim_out(end, :)];
    timing = [timing; output.info.cpuTime];
    status = [status; output.info.status];
    num_iters = [num_iters; output.info.QP_iter];
    
end

system('rm acado_MPCstep.mexmaci64')

end

