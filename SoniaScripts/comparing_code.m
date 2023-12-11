% Setup

% par_og = Richardson2000FullAxon();
par_crsp = parameters_Richardson();
par_me = Cullen2018CortexAxonJPN_MTR();
par_me1 = parameters_corticalaxon();

par_me2 = parameters_Richardson_myel();
% par_me2.channels(4).cond.value = 0 * par_me2.channels(4).cond.value; % toggle
% on/off to see influence of channels under myelin

%% test - same as run_Richardson_default
par = par_me2;
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
title(['Conduction velocity (m/s): ' num2str(CV)])
xlabel('Time')
ylabel('Voltage')


%% to run Cullen model
% par = Cullen2018CortexAxon(); % no active channels added in the internode region
% 
% % Run sham simulations.
% [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model(par);
% 
% 
% cv_between_nodes = [6, 13];
% internode_length = unitsabs(par.intn.geo.length.units) * par.intn.geo.length.value.ref + ...
%     unitsabs(par.node.geo.length.units) * par.node.geo.length.value.ref;
% 
% [~, tmax2] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(2)));
% [~, tmax1] = max(MEMBRANE_POTENTIAL(:, cv_between_nodes(1)));
% 
% CV = diff(cv_between_nodes) * internode_length / ((tmax2 - tmax1) * dt);
% fprintf('Conduction velocity: %.8f m/s\n', CV);
% figure
% hold on
% plot(TIME_VECTOR,MEMBRANE_POTENTIAL);
% title(['Conduction velocity (m/s): ' num2str(CV)])
% xlabel('time')
% ylabel('voltage')
% hold off
