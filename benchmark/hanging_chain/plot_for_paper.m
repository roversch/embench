
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

markerstyles = {'.', 'x', 'd', 'v', 'o'};

figure(2); clf;
for j=1:numel(solvers)
    distance_to_ref = max(abs(X.(solvers{j}) - repmat(x_ref.', num_sim_iters+1, 1)).');
    semilogy(time_vector, distance_to_ref, 'LineWidth', 2, 'Marker', markerstyles{j}, 'LineStyle', '-.');
    hold on
end
legend(solvers, 'Location', 'SouthWest');
set(gca,'FontSize',30);
set(findall(gcf,'type','text'),'FontSize',30);
xlabel('Simulation time $[\mathrm{s}]$', 'Interpreter', 'Latex')
ylabel('$\|\overline{x}_0 - x_\mathrm{ref}\|$', 'Interpreter', 'Latex')

figure(3); clf;
for j=1:numel(solvers)
   semilogy(time_vector(1:end-1), timing.(solvers{j}), 'LineWidth', 2, 'Marker', markerstyles{j}, 'LineStyle', '-.');
   hold on;
end
legend(solvers, 'Location', 'SouthWest');

set(gca,'FontSize',30);
set(findall(gcf,'type','text'),'FontSize',30);


