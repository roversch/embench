
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);

figure(1); clf;
for i=1:4
    subplot(4, 2, i); hold on;
    for j=1:numel(solvers)
        plot(time_vector, X.(solvers{j})(:, i))
    end
    title(['x ' num2str(i)])
end

subplot(4, 2, 5); hold on;
for j=1:numel(solvers)
    stairs(time_vector(1:end-1), U.(solvers{j}))
    title(['u'])
end

subplot(4, 2, 6); hold on;
for j=1:numel(solvers)
    stairs(time_vector(1:end-1), status.(solvers{j}))
    title('status')
end

subplot(4, 2, 7); hold on;
for j=1:numel(solvers)
   semilogy(timing.(solvers{j}));
   title('cpu time')
end
legend(solvers);

subplot(4, 2, 8); hold on;
for j=1:numel(solvers)
   stairs(num_iters.(solvers{j}));
   title('qp iters')
end

fprintf([repmat('-', 1, 80), '\n']);
fprintf(['(ms)\t\tmean\t\t\tmin\t\t\tmax\t\t', '\n']);
fprintf([repmat('-', 1, 80), '\n']);
for j=1:numel(solvers)
   fprintf(['%s\t\t', '%.2f', '\t\t\t', '%.2f', '\t\t\t', '%.2f', '\t\t\n'], ...
       solvers{j}, 1000*mean(timing.(solvers{j})), ...
       1000*min(timing.(solvers{j})), 1000*max(timing.(solvers{j})));
end
