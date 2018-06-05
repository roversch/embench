
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

markerstyles = {'.', 'x', 'd', 'v', 'o', 's'};

figure(2); clf;
for j=1:numel(solvers)
    
    XU = [X.(solvers{j}), [U.(solvers{j}); zeros(1, 3)]];
    
    XUerr = XU - repmat([x_ref.', zeros(1, 3)], num_sim_iters+1, 1);
    
    distance_to_ref = diag(XUerr * W * XUerr.');
    semilogy(time_vector, distance_to_ref, 'LineWidth', 2, 'Marker', markerstyles{j}, 'LineStyle', '-.');
    hold on
end
legend(solvers);
set(gca,'FontSize',30);
set(findall(gcf,'type','text'),'FontSize',30);
xlabel('Simulation time $[\mathrm{s}]$', 'Interpreter', 'Latex')
ylabel('$\|x_0-x_\mathrm{ref}\|_Q^2 + \|u_0\|_R^2$', 'Interpreter', 'Latex')

figure(3); clf;
for j=1:numel(solvers)
   semilogy(time_vector(1:end-1), timing.(solvers{j}), 'LineWidth', 2, 'Marker', markerstyles{j}, 'LineStyle', '-.');
   hold on;
end
legend(solvers);

set(gca,'FontSize',30);
set(findall(gcf,'type','text'),'FontSize',30);

xlabel('Simulation time $[\mathrm{s}]$', 'Interpreter', 'Latex')
ylabel('CPU time $[\mathrm{s}]$', 'Interpreter', 'Latex')

