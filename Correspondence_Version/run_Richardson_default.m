% example for running the model with the default Richardson et al. (2000)
% parameters
par = parameters_Richardson;
V = model(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
T = (0:dt:dt*(size(V, 3)-1));
node_idx = [par.geo.nodeSegments{:}];

% voltage in nodes
Vnodes = reshape(V(node_idx, 1, :), length(node_idx), [])';

% plot voltage in nodes
figure, plot(T, Vnodes)

% conduction velocity
cv_between_nodes = [6, 13];
internode_length = unitsabs(par.intn.geo.length.units) * par.intn.geo.length.value.ref + ...
    unitsabs(par.node.geo.length.units) * par.node.geo.length.value.ref;

[~, tmax2] = max(Vnodes(:, cv_between_nodes(2)));
[~, tmax1] = max(Vnodes(:, cv_between_nodes(1)));

CV = diff(cv_between_nodes) * internode_length / ((tmax2 - tmax1) * dt);
fprintf('Conduction velocity: %.8f m/s\n', CV);