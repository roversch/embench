function [states, controls, timing, status, num_iters] = grampc_run(num_free_masses, N, Ts, W, WN, umax, x_ref, u_0, num_sim_iters, integrator_fun)

if (num_free_masses ~= 1)
    nx = (2*num_free_masses+1)*3;
    nu = 3;
    states = zeros(num_sim_iters, nx);
    controls = zeros(num_sim_iters, nu);
    timing = zeros(num_sim_iters, 1);
    status = ones(num_sim_iters, 1);
    num_iters = zeros(num_sim_iters, 1);
    return
end

probfct = 'grampc_hanging_chain_model.c';

copyfile(['../model/', probfct], probfct);

grampc_root_path = '/Users/robin/GRAMPC/';

addpath([grampc_root_path 'matlab/mfiles']);

% grampc_make_toolbox(grampc_root_path);
% grampc_make_probfct(grampc_root_path, probfct);

% Parameter definition
% Initial values and setpoints of the states
user.param.x0    = x_ref;
user.param.xdes  = x_ref;

% Initial values, setpoints and limits of the inputs
user.param.u0    = zeros(3, 1);
user.param.udes  = zeros(3, 1);
user.param.umax  = +umax;
user.param.umin  = -umax;

% Time variables
user.param.Thor  = N*Ts;           % Prediction horizon

user.param.dt    = Ts;          % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

% Option definition
% Basic algorithmic options
user.opt.Nhor        = N;      % Number of steps for the system integration
user.opt.TerminalCost = 'on';
user.opt.ShiftControl = 'off';
user.opt.Integrator = 'ruku45';
user.opt.MaxMultIter = 2;
user.opt.MaxGradIter = 23;
user.opt.EqualityConstraints = 'off';
user.opt.InequalityConstraints = 'off';
user.opt.TerminalEqualityConstraints = 'off';
user.opt.TerminalInequalityConstraints = 'off';
user.opt.ConvergenceCheck = 'on';

% Constraints tolerances
% user.opt.ConstraintsAbsTol = 1e-3*[1e-1 1 1];

% optional settings for a better performance
% user.opt.LineSearchMax = 1e1;
% user.opt.LineSearchInit = 1e-1;
% user.opt.LineSearchExpAutoFallback = 'off';
% user.opt.PenaltyMin = 4e1; % Comment line 76 (grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);) to use this option value

% Userparameter definition 
% e.g. system parameters or weights for the cost function
userparam = [diag(W); diag(WN)];

% Grampc initialization
grampc = CmexFiles.grampc_init_Cmex(userparam);

% Update grampc struct while ensuring correct data types
grampc = grampc_update_struct_grampc(grampc,user);

% Estimate and set PenaltyMin
grampc = CmexFiles.grampc_estim_penmin_Cmex(grampc,1);

Tsim = num_sim_iters * Ts;

CmexFiles.grampc_printopt_Cmex(grampc);
CmexFiles.grampc_printparam_Cmex(grampc);

% init solution structure
vec = grampc_init_struct_sol(grampc, Tsim);

states = x_ref.';
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
        
    if i <= 5
        controls = [controls; u_0.'];
    else
        controls = [controls; vec.u(:, i).'];
    end
    
    [~, sim_out] = integrator_fun(states(end, :).', controls(end, :).');
    
    states = [states; sim_out(end, :)];
    timing = [timing; vec.CPUtime(i)/1000];
    status = [status; grampc.sol.status];
    num_iters = [num_iters; vec.iter(1, i)];

    % update iteration counter
    vec.x(:,i+1) = states(end, :).';
    i = i + 1;
    
end

rmpath([grampc_root_path 'matlab/mfiles']);

end

