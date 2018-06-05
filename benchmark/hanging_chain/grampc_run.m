function [states, controls, timing, status, num_iters] = grampc_run(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

rng(0);

probfct = 'probfct_NLCHAIN_4.c';

copyfile(['../model/', probfct], probfct);

grampc_root_path = '/Users/robin/GRAMPC/';

addpath([grampc_root_path 'matlab/mfiles']);

grampc_make_toolbox(grampc_root_path);
grampc_make_probfct(grampc_root_path, probfct);

%% Parameter definition
% Initial values and setpoints of the states

xdim = 3*(2*num_free_masses+1);
user.param.x0    = x_0.';
user.param.xdes  = zeros(1,xdim);

% Initial values, setpoints and limits of the inputs
user.param.u0    = [0.0,0.0,0.0];
user.param.udes  = [0.0,0.0,0.0];
user.param.umax  = +umax;
user.param.umin  = -umax;

% Time variables
user.param.Thor  = N*Ts;         % Prediction horizon
user.param.Nhor  = N;

user.param.dt    = Ts;         % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

%% Option definition
% Basic algorithmic options
user.opt.MaxMultIter = 5;           % Maximum number of augmented Lagrange iterations
% user.opt.MaxGradIter = 20;           % Maximum number of gradient iterations
user.opt.ShiftControl = 'off';
user.opt.ConvergenceCheck = 'on';

% Cost integration
user.opt.TerminalCost = 'on';

% System integration
user.opt.Integrator = 'ruku45';
user.opt.IntegratorRelTol = 1e-3;
user.opt.IntegratorAbsTol = 1e-4;

% optional settings for a better performance
% user.opt.LineSearchType ='adaptive';
% user.opt.LineSearchMax = 1e2;
% user.opt.LineSearchInit = 1e-1;
% user.opt.LineSearchExpAutoFallback = 'off';

%% Userparameter definition 
% e.g. system parameters or weights for the cost function
pCost = [25, 2.5, 0.1, 10];
userparam = pCost;

% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

Tsim = num_sim_iters * Ts;

CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

avg_timing = 0;

for irepeat=1:nrepeat

% init solution structure
vec = grampc_init_struct_sol(grampc, Tsim);

states = x_0.';
controls = [];
timing = [];
status = [];
num_iters = [];

i = 1;
while 1
    % set current time and current state
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'t0',vec.t(i));
    grampc = CmexFiles.grampc_setparam_Cmex(grampc,'x0',vec.x(:,i));
    
    % run MPC and save results
    [grampc,vec.CPUtime(i)] = CmexFiles.grampc_run_Cmex(grampc);
    vec = grampc_update_struct_sol(grampc, vec, i);
    
    % print solver status
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status, 'Debug');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
        
    % check for end of simulation
    if i+1 > length(vec.t)
        break;
    end
    
    % simulate system
        
    if i >= num_sim_iters/2 && i < num_sim_iters/2 + 5
        controls = [controls; u_0.'];
    else
        controls = [controls; vec.u(:, i).'];
    end
        
    [~, sim_out] = integrator_fun(states(end, :).', controls(end, :).');
    
    states = [states; sim_out(end, :) + sigma*randn(1, xdim)];
    timing = [timing; vec.CPUtime(i)/1000];
    status = [status; grampc.sol.status];
    num_iters = [num_iters; vec.iter(1, i)];
    
    % update iteration counter
    vec.x(:,i+1) = states(end, :).';
    i = i + 1;
    
end

avg_timing = avg_timing + timing;

end

timing = avg_timing / nrepeat;

rmpath([grampc_root_path 'matlab/mfiles']);

end

