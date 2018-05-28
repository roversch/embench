function [states, controls, timing, status, num_iters] = grampc_run(N, Ts, W, WN, Fmax, x0, num_sim_iters)

probfct = 'grampc_pendulum_model.c';

copyfile(['../model/', probfct], probfct);

grampc_root_path = '/Users/robin/GRAMPC/';

addpath([grampc_root_path 'matlab/mfiles']);

grampc_make_toolbox(grampc_root_path);
grampc_make_probfct(grampc_root_path, probfct);

% Parameter definition
% Initial values and setpoints of the states
user.param.x0    = x0;
user.param.xdes  = [0; 0; 0; 0];

% Initial values, setpoints and limits of the inputs
user.param.u0    = 0.0;
user.param.udes  = 0.0;
user.param.umax  = +Fmax;
user.param.umin  = -Fmax;

% Time variables
user.param.Thor  = N*Ts;           % Prediction horizon

user.param.dt    = Ts;          % Sampling time
user.param.t0    = 0.0;         % time at the current sampling step

% Option definition
% Basic algorithmic options
user.opt.Nhor        = N;      % Number of steps for the system integration
user.opt.TerminalCost = 'on';
user.opt.ShiftControl = 'off';
user.opt.MaxMultIter = 20;
user.opt.MaxGradIter = 20;

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

states = x0.';
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
    printed = CmexFiles.grampc_printstatus_Cmex(grampc.sol.status,'Error');
    if printed
        fprintf('at simulation time %f.\n --------\n', vec.t(i));    
    end
    
    % check for end of simulation
    if i+1 > length(vec.t)
        break;
    end
    
    % simulate system
    vec.x(:,i+1) = grampc.sol.xnext;
        
    states = [states; vec.x(:, i+1).'];
    controls = [controls; vec.u(:, i).'];
    timing = [timing; vec.CPUtime(i)/1000];
    status = [status; grampc.sol.status];
    num_iters = [num_iters; vec.iter(1, i)];

    % update iteration counter
    i = i + 1;
    
end

rmpath([grampc_root_path 'matlab/mfiles']);

end

