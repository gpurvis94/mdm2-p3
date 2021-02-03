clc;
close all;


% initialise variables
a = 0.2;                    % disease transmission probability
b = 2;                      % contact rate
beta = a*b;                 % infection rate
infectivePeriod = 7;        % how long infection lasts
gamma = 1/infectivePeriod;  % recovery rate

t0 = 0;         % start time
tf = 56;        % end time
h = 1;          % time step

N = 60;     %population size
I0 = 1;        %infected population
R = 0;        %recovered population
S = N - I0;    %susceptible population

%initial conditions for ODE
y0 = [S, I0, R]; 



% Solving ODE
[t,y] = ode45(@(t,y) SIRRHS(t,y,N,beta,gamma), [t0 tf], y0);

%plotting solution from ode45
plot(t, y(:,1), '-g', 'LineWidth', 2)
hold on
plot(t, y(:,2), '-r', 'LineWidth', 2)
plot(t, y(:,3), '-b', 'LineWidth', 2)
legend('S','I','R');
title('SIR Model', 'FontSize', 20);
xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
hold off

%generating networks
smallWorld = full(small_world(N, 8, 0.001));
scaleFree = full(scale_free(N, 6, 2));
randomGraph = full(random_graph(N, 0.2));
completeGraph = full(complete_graph(N));

%calculating clusters
smallWorldCluster = global_clustering_coefficient(smallWorld);
scaleFreeCluster = global_clustering_coefficient(scaleFree);
randomCluster = global_clustering_coefficient(randomGraph);
completeCluster = global_clustering_coefficient(completeGraph);

%running the model
[SVec1, IVec1, RVec1, smallWorldData] = networkSIR(N, a, b, infectivePeriod, smallWorld, tf, I0);
[SVec2, IVec2, RVec2, scaleFreeData] = networkSIR(N, a, b, infectivePeriod, scaleFree, tf, I0);
[SVec3, IVec3, RVec3, randomData] = networkSIR(N, a, b, infectivePeriod, randomGraph, tf, I0);
[SVec4, IVec4, RVec4, completeData] = networkSIR(N, a, b, infectivePeriod, completeGraph, tf, I0);
%Data Matrices
%row 1 is the node number
%row 2 is infection status
%row 3 is the day of infection
%row 4 is the local cluster coefficient
%row 5 is the degree of each node
%row 6 is the shortest distance from node to nearest patient0
%row 7 is the number of nodes infected by this node

%%%%%%%%%
%Analyse existing networks
%[A, rows] = mmread('socfb-Amherst41.mtx');
%[SVec5, IVec5, RVec5, ameherstData] = networkSIR(rows, a, b, infectivePeriod, A, tf, I0);
%amherstCluster = global_clustering_coefficient(A);

%figure;
%tiledlayout(2,1);
%nexttile;
%plot((0:length(SVec5)-1), SVec5, '-g', 'LineWidth', 2)
%hold on
%plot((0:length(IVec5)-1), IVec5, '-r', 'LineWidth', 2)
%plot((0:length(RVec5)-1), RVec5, '-b', 'LineWidth', 2)
%legend('S','I','R');
%title(sprintf('SIR Network Model (Amherst FB) GCC = %0.5f',amherstCluster), 'FontSize', 20);
%xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
%grid on;
%ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
%nexttile;
%G = graph(A);
%p=plot(graph(A))
%title('Amherst FB network');
%G.Nodes.NodeColors = degree(G);
%p.NodeCData = G.Nodes.NodeColors;
%colorbar;

%%%%%%%%


%plotting solutions for networked model
figure;
tiledlayout(2,2);
nexttile;
plot((0:length(SVec1)-1), SVec1, '-g', 'LineWidth', 2)
hold on
plot((0:length(IVec1)-1), IVec1, '-r', 'LineWidth', 2)
plot((0:length(RVec1)-1), RVec1, '-b', 'LineWidth', 2)
legend('S','I','R');
title(sprintf('SIR Network Model (small-world) GCC = %0.5f',smallWorldCluster), 'FontSize', 20);
xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
hold off
nexttile;
plot((0:length(SVec2)-1), SVec2, '-g', 'LineWidth', 2)
hold on
plot((0:length(IVec2)-1), IVec2, '-r', 'LineWidth', 2)
plot((0:length(RVec2)-1), RVec2, '-b', 'LineWidth', 2)
legend('S','I','R');
title(sprintf('SIR Network Model (scale-free) GCC = %0.5f',scaleFreeCluster), 'FontSize', 20);
xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
hold off
nexttile;
plot((0:length(SVec3)-1), SVec3, '-g', 'LineWidth', 2)
hold on
plot((0:length(IVec3)-1), IVec3, '-r', 'LineWidth', 2)
plot((0:length(RVec3)-1), RVec3, '-b', 'LineWidth', 2)
legend('S','I','R');
title(sprintf('SIR Network Model (random) GCC = %0.5f',randomCluster), 'FontSize', 20);
xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
hold off
nexttile;
plot((0:length(SVec4)-1), SVec4, '-g', 'LineWidth', 2)
hold on
plot((0:length(IVec4)-1), IVec4, '-r', 'LineWidth', 2)
plot((0:length(RVec4)-1), RVec4, '-b', 'LineWidth', 2)
legend('S','I','R');
title(sprintf('SIR Network Model (complete) GCC = %0.5f',completeCluster), 'FontSize', 20);
xlabel('time', 'FontSize', 20),ylabel('people', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
hold off

%show network graphs
figure;
tiledlayout(2,2);
nexttile;
G = graph(smallWorld);
p=plot(graph(smallWorld))
title('small-world network');
G.Nodes.NodeColors = degree(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar;
nexttile;
G = graph(scaleFree);
p=plot(graph(scaleFree))
title('scale-free network');
G.Nodes.NodeColors = degree(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar;
nexttile;
G = graph(scaleFree);
p=plot(graph(randomGraph))
title('random network');
G.Nodes.NodeColors = degree(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar;
nexttile;
G = graph(completeGraph);
p=plot(graph(completeGraph))
title('complete network');
G.Nodes.NodeColors = degree(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar;

%day node got infected against distance to P0
figure;
tiledlayout(2,2);
nexttile;
scatter(smallWorldData(6,:),smallWorldData(3,:));
title('Time of infection against steps to P0 (small-world)', 'FontSize', 20);
xlabel('distance to P0', 'FontSize', 20),ylabel('time of infection', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(scaleFreeData(6,:),scaleFreeData(3,:));
title('Time of infection against steps to P0 (scale-free)', 'FontSize', 20);
xlabel('distance to P0', 'FontSize', 20),ylabel('time of infection', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(randomData(6,:),randomData(3,:));
title('Time of infection against steps to P0 (random)', 'FontSize', 20);
xlabel('distance to P0', 'FontSize', 20),ylabel('time of infection', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(completeData(6,:),completeData(3,:));
title('Time of infection against steps to P0 (complete)', 'FontSize', 20);
xlabel('distance to P0', 'FontSize', 20),ylabel('time of infection', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;

%No of infections caused against LCC
figure;
tiledlayout(2,2);
nexttile;
scatter(smallWorldData(4,:),smallWorldData(7,:));
title('No of infections caused against LCC (small-world)', 'FontSize', 20);
xlabel('LCC', 'FontSize', 20),ylabel('No of infections caused', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(scaleFreeData(4,:),scaleFreeData(7,:));
title('No of infections caused against LCC (scale-free)', 'FontSize', 20);
xlabel('LCC', 'FontSize', 20),ylabel('No of infections caused', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(randomData(4,:),randomData(7,:));
title('No of infections caused against LCC (random)', 'FontSize', 20);
xlabel('LCC', 'FontSize', 20),ylabel('No of infections caused', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(completeData(4,:),completeData(7,:));
title('No of infections caused against LCC (complete)', 'FontSize', 20);
xlabel('LCC', 'FontSize', 20),ylabel('No of infections caused', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;

%number of infections caused against node degree
figure;
tiledlayout(2,2);
nexttile;
scatter(smallWorldData(5,:),smallWorldData(7,:));
title('Number of infections caused against degree (small-world)', 'FontSize', 20);
xlabel('degree', 'FontSize', 20),ylabel('Number of infections', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(scaleFreeData(5,:),scaleFreeData(7,:));
title('Number of infections caused against degree (scale-free)', 'FontSize', 20);
xlabel('degree', 'FontSize', 20),ylabel('Number of infections', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(randomData(5,:),randomData(7,:));
title('Number of infections caused against degree (random)', 'FontSize', 20);
xlabel('degree', 'FontSize', 20),ylabel('Number of infections', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
nexttile;
scatter(completeData(5,:),completeData(7,:));
title('Number of infections caused against degree (complete)', 'FontSize', 20);
xlabel('degree', 'FontSize', 20),ylabel('Number of infections', 'FontSize', 20);
grid on;
ax = gca; ax.YAxis.FontSize = 15; ax.XAxis.FontSize = 15;
