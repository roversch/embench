
time_vector = linspace(0, num_sim_iters*Ts, num_sim_iters+1);

solvers = fieldnames(X);
suboptimality = struct();

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

markerstyles = {'.', 'x', 'd', 'v', 'o', 's'};

figure(2); clf; handles = [];
for j=1:numel(solvers)
    
    XU = [X.(solvers{j}), [U.(solvers{j}); zeros(1, 3)]];
    
    XUerr = XU - repmat([x_ref.', zeros(1, 3)], num_sim_iters+1, 1);
    
    distance_to_ref = diag(XUerr * W * XUerr.');
    if j == 1
        % ipopt
        ipopt_suboptimality = distance_to_ref;
    end
    suboptimality.(solvers{j}) = abs((cumsum(distance_to_ref)-cumsum(ipopt_suboptimality))./cumsum(ipopt_suboptimality));
    handles = [handles, semilogy(time_vector, suboptimality.(solvers{j}), 'LineWidth', 2, 'Marker', markerstyles{j}, 'LineStyle', '-.')];
    hold on;
end
legend(handles(2:end), solvers{2:end});
set(gca,'FontSize',30);
set(findall(gcf,'type','text'),'FontSize',30);
xlabel('Simulation time $[\mathrm{s}]$', 'Interpreter', 'Latex')
ylabel('relative cumulative suboptimality', 'Interpreter', 'Latex')

grid on

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
grid on


figure(4); clf; 
for j=1:numel(solvers)
   subopt = suboptimality.(solvers{j});
   subopt = subopt(end);
   loglog(max(timing.(solvers{j})), subopt, 'ok', 'MarkerSize', 10, 'LineWidth', 4);
   hold on;
end

xlim([1e-3, 2e-2])
ylim([1e-5, 1])
grid on;

xlabel('max CPU time $[\mathrm{s}]$', 'Interpreter', 'Latex')
ylabel('relative cumulative suboptimality', 'Interpreter', 'Latex')

figure1 = gcf;

% Create textbox
annotation(figure1,'textbox',...
    [0.644947866258475 0.578034050358245 0.0915110132158591 0.0508130081300824],...
    'String','\texttt{viatoc}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.565337431776313 0.860351874599121 0.0915110132158589 0.0508130081300824],...
    'String','\texttt{falcopt}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.21083125121783 0.294332284928956 0.091511013215859 0.0508130081300822],...
    'String','\texttt{acados}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.490841044091028 0.300748356591569 0.091511013215859 0.0508130081300822],...
    'String','\texttt{acado}',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
    [0.209585198499028 0.751995546335499 0.0915110132158588 0.0508130081300824],...
    'String',{'\texttt{grampc}'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',40,...
    'FitBoxToText','off');

set(gca,'FontSize',40);
set(findall(gcf,'type','text'),'FontSize',30);
