function [states, controls, timings, status] = acado_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

clear GLOBAL

DifferentialState p theta v dtheta;
Control Force;

% Export of an optimization routine:
acadoSet('problemname', 'pendulum');

ocp = acado.OCP(0.0, N*Ts, N);

h = [p; theta; v; dtheta; Force];
hN = [p; theta; v; dtheta];

ocp.minimizeLSQ(W, h);            % stage cost
ocp.minimizeLSQEndTerm(WN, hN);   % terminal cost

Fmin = -Fmax; 
Fmax = +Fmax;

ode = pendulum_ode([p; theta; v; dtheta], Force);
ocp.setModel(ode);  % pass the ODE model

ocp.subjectTo(Fmin <= Force <= Fmax);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',        N                   );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);

mpc.exportCode('pendulum_export');

global ACADO_;
copyfile([ACADO_.pwd '/../../external_packages/qpoases'], 'pendulum_export/qpoases')

cd pendulum_export
make_acado_solver('../acado_MPCstep')
cd ..

% Parameters
input.x = repmat(x0.', N+1, 1);
input.u = zeros(N, 1);

input.y = zeros(N, 5);
input.yN = zeros(1, 4);

input.lbValues = -Fmax*ones(N, 1);
input.ubValues = +Fmax*ones(N, 1);

input.shifting.strategy = 0;

states = x0.';
controls = [];
timings = [];
status = [];

for i=1:num_sim_iters
    
    % Solve NMPC OCP
    input.x0 = states(end, :);
    
    
    kkt = inf;
    running_time = 0;
    while (kkt > 1e-10)
        output = acado_MPCstep(input);
        input.x = output.x;
        input.u = output.u;
        kkt = output.info.kktValue;
        running_time = running_time + output.info.cpuTime;
    end
       
    % Save the MPC step
    states = [states; output.x(2, :)];
    controls = [controls; output.u(1, :)];
    timings = [timings; running_time];
    status = [status; output.info.status];
    
    input.x = output.x;
    input.u = output.u;
    
end  % while

end  % function

