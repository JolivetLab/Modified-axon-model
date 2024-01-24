% Setup
clear
% par_og = Richardson2000FullAxon();
par_crsp = parameters_Richardson();
par_default = Cullen2018CortexAxon();
par_me = Cullen2018CortexAxonJPN_MTR();
par_me1 = parameters_corticalaxon();
par_me2 = parameters_Richardson_myel();
% on/off to see influence of channels under myelin

%% test - same as run_Richardson_default
par = parameters_Richardson_myel();
V = model(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
T = (0:dt:dt*(size(V, 3)-1));
node_idx = [par.geo.nodeSegments{:}];

% voltage in nodes
Vnodes = reshape(V(node_idx, 1, :), length(node_idx), [])';

% conduction velocity
cv_between_nodes = [6, 13];
internode_length = unitsabs(par.intn.geo.length.units) * par.intn.geo.length.value.ref + ...
    unitsabs(par.node.geo.length.units) * par.node.geo.length.value.ref;

[~, tmax2] = max(Vnodes(:, cv_between_nodes(2)));
[~, tmax1] = max(Vnodes(:, cv_between_nodes(1)));

CV = diff(cv_between_nodes) * internode_length / ((tmax2 - tmax1) * dt);
fprintf('Conduction velocity: %.8f m/s\n', CV);

% plot voltage in nodes
figure
plot(T, Vnodes)
title(['Lee - Conduction velocity (m/s): ' num2str(CV)])
xlabel('Time')
ylabel('Voltage')


%% to run modified model
par = Cullen2018CortexAxonJPN_MTR();

% Run sham simulations.
[MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
% T = (0:dt:dt*(size(MEMBRANE_POTENTIAL, 3)-1));

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

%% to run Cullen model
par = Cullen2018CortexAxon();

% Run sham simulations.
[MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
T = (0:dt:dt*(size(MEMBRANE_POTENTIAL, 3)-1));
% node_idx = [par.geo.nodeSegments{:}];

cv_between_nodes = [6, 13];

[~, tmax2] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(2)));
[~, tmax1] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(1)));

CV = diff(cv_between_nodes) * INTERNODE_LENGTH(1) / ((tmax2 - tmax1) * dt);
fprintf('Conduction velocity: %.8f m/s\n', CV);
figure
hold on
plot(TIME_VECTOR,MEMBRANE_POTENTIAL);
title(['Cullen - Conduction velocity (m/s): ' num2str(CV)])
xlabel('time')
ylabel('voltage')
hold off
