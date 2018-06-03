function [states, controls, timing, status, num_iter] = run_viatoc(num_free_masses, N, Ts, W, WN, umax, x_ref, u_0, num_sim_iters, integrator_fun)

x_ref_traj = repmat(x_ref, 1, N);

sim('run_viatoc_sim');

states = states.data;
controls = controls.data(1:end-1, :);
timing = zeros(size(controls, 1), 1);
status = zeros(size(controls, 1), 1);
num_iter = zeros(size(controls, 1), 1);

end

