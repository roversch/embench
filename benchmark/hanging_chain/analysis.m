
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);

distance_to_ref = struct();

figure(1); clf;
subplot(4, 2, 1);
for j=1:numel(solvers)
    
    XU = [X.(solvers{j}), [U.(solvers{j}); zeros(1, 3)]];
    
    XUerr = XU - repmat([x_ref.', zeros(1, 3)], num_sim_iters+1, 1);
    
    distance_to_ref.(solvers{j}) = diag(XUerr * W * XUerr.');
    title('weighted distance to reference');
    semilogy(time_vector, distance_to_ref.(solvers{j}));
    hold on
end

subplot(4, 2, 2);
for j=1:numel(solvers)
    title('distance to IPOPT');
    distance_to_ipopt = max(abs(U.ipopt - U.(solvers{j})), [], 2);
    stairs(time_vector(1:end-1), distance_to_ipopt);
    hold on
end
plot(time_vector(1:end-1), 0*distance_to_ipopt + 0.2, 'k--', 'LineWidth', 2);

for i=1:3
    subplot(4, 2, 2+i); hold on;
    title(['u' num2str(i)]);
    for j=1:numel(solvers)
        stairs(time_vector(1:end-1), U.(solvers{j})(:, i))
    end
end

subplot(4, 2, 6); hold on;
for j=1:numel(solvers)
    title('solver status');
    stairs(time_vector(1:end-1), status.(solvers{j}))
end

subplot(4, 2, 7);
for j=1:numel(solvers)
   title('solver timing');
   semilogy(timing.(solvers{j}));
   hold on;
end
legend(solvers);

subplot(4, 2, 8);
for j=1:numel(solvers)
   title('solver iterations');
   semilogy(num_iters.(solvers{j}));
   hold on;
end

fprintf([repmat('-', 1, 100), '\n']);
fprintf(['(ms)\t\tmedian\t\t\tmin\t\t\tmax\t\t\tcost', '\n']);
fprintf([repmat('-', 1, 100), '\n']);
for j=1:numel(solvers)
   fprintf(['%s\t\t', '%.2f', '\t\t\t', '%.2f', '\t\t\t', '%.2f', '\t\t\t', '%.2f', '\n'], ...
       solvers{j}, 1000*median(timing.(solvers{j})), ...
       1000*min(timing.(solvers{j})), 1000*max(timing.(solvers{j})), ...
       max(max(abs(((U.(solvers{j}) - U.ipopt))))));
end

figure(2);
for j=1:numel(solvers)
   title('solver timing');
   semilogy(timing.(solvers{j}));
   hold on;
end
legend(solvers);
print(gcf, '-dpng', '-r100', 'CPU_timing');

