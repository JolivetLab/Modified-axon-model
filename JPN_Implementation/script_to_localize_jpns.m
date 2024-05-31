clear

%% plot coordinates
par = Cullen2018CortexAxonJPNlocalized_MTR();

nodes = cell2mat(par.geo.nodeSegments);
intn = cell2mat(par.geo.internodeSegments);
jpns = cell2mat(par.geo.juxtaparanodeSegments);

x_nodes = zeros(size(par.geo.nodeSegments));
x_intn = zeros(size(par.geo.internodeSegments));
x_jpn = zeros(size(par.geo.juxtaparanodeSegments));

figure;
hold on
plot(nodes, x_nodes, 'm.', 'LineWidth', 1, 'MarkerSize', 15);
plot(intn, x_intn, 'b.', 'LineWidth', 1, 'MarkerSize', 15);
plot(jpns, x_jpn, 'g.', 'LineWidth', 1, 'MarkerSize', 15);
title(['magenta - node, blue - internode, green - juxtaparanode']);
xlim([1 160]);
ylim([-0.5 0.5])
xlabel('segment');
% title('Overlap of Coordinates');
grid on;
%% RUN SIMULATION
par = Cullen2018CortexAxonJPNlocalized_MTR();
par.sim.dt.value = 1;
[MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;

cv_between_nodes = [6, 13];

internode_length = unitsabs(par.intn.geo.length.units) * par.intn.geo.length.value.ref + ...
    unitsabs(par.node.geo.length.units) * par.node.geo.length.value.ref;

[~, tmax2] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(2)));
[~, tmax1] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(1)));

CV = diff(cv_between_nodes) * internode_length / ((tmax2 - tmax1) * dt);
fprintf('Conduction velocity: %.8f m/s\n', CV);
figure
hold on
plot(TIME_VECTOR,MEMBRANE_POTENTIAL);
title(['My version - Conduction velocity (m/s): ' num2str(CV)])
xlabel('time')
ylabel('voltage')
hold off
%% 
% par = Cullen2018CortexAxonJPNlocalized_MTR;
% CalculateLeakConductanceJPN_MTR(par);