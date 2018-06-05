function [states, controls, timing, status, num_iter] = run_viatoc(num_free_masses, N, Ts, W, WN, umax, x_ref, x_0, u_0, num_sim_iters, integrator_fun, sigma, nrepeat)

timing = 0;
for i=1:nrepeat
    sim('viatoc_run_sim');
    timing = timing + viatoc_timing.data(1:end-1,:);
end

states = states.data;
controls = controls.data(1:end-1, :);
timing = timing / nrepeat;
status = zeros(size(controls, 1), 1);
num_iter = zeros(size(controls, 1), 1);

end

