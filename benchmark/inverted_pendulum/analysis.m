
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);

figure(1); clf;
for i=1:4
    subplot(3, 2, i); hold on;
    for j=1:numel(solvers)
        plot(time_vector, X.(solvers{j})(:, i))
    end
end

subplot(3, 2, 5); hold on;
for j=1:numel(solvers)
    stairs(time_vector(1:end-1), U.(solvers{j}))
end

subplot(3, 2, 6); hold on;
for j=1:numel(solvers)
    stairs(time_vector(1:end-1), status.(solvers{j}))
end
