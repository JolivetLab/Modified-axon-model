%% TEST SIMULATIONS
activeChannel   = JPN_SlowKv11; save('D:/0_UNI/PRO2 Myelin/Model/SavedParameters/ActiveChannels/JPN_SlowKv11.mat','activeChannel');

par = Cullen2018CortexAxon_jpn_test();
% par.intn.elec.act.erev.value = -10;
% par.intn.elec.act.cond.value.ref = 0.0;
% par.intn.elec.act.cond.value.vec = zeros(4,4);
% set input to 0 - should see flatlines or potentially a small drift
par_ctr = Cullen2018CortexAxon(); % no active channels added in the internode region


% Set temperature.
par.sim.temp = 37;
par_ctr.sim.temp = 37;
    
% Adjust simulation time.
par.sim.tmax.value = 15;
par_ctr.sim.tmax.value = 15;
% par.sim.dt.value = 0.005;
    
% Run sham simulations.
% [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = Model_jpn(par);
[MEMBRANE_POTENTIAL_ctr, INTERNODE_LENGTH_ctr, TIME_VECTOR_ctr] = Model(par_ctr);
[MEMBRANE_POTENTIAL_ctr2, INTERNODE_LENGTH_ctr2, TIME_VECTOR_ctr2] = Model(par); % Model_copy_nojpn(par);

figure
hold on
plot(TIME_VECTOR_ctr,mean(MEMBRANE_POTENTIAL_ctr,2),'-k','LineWidth',1.5, 'DisplayName','Modified Axon');
% plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-b', 'DisplayName','Modified Model');
plot(TIME_VECTOR_ctr2,mean(MEMBRANE_POTENTIAL_ctr2,2),'--r','LineWidth',1.5, 'DisplayName', 'Base Axon')
legend
xlabel('time')
ylabel('voltage')
hold off

figure
hold on
plot(TIME_VECTOR_ctr,MEMBRANE_POTENTIAL_ctr,'-k');
% plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-b', 'DisplayName','Modified Model');
plot(TIME_VECTOR_ctr2,MEMBRANE_POTENTIAL_ctr2,'--r')
xlabel('time')
ylabel('voltage')
hold off

%% GEOMETRY FIDDLING
% nis = par.geo.nintseg;
% nns = par.geo.nnodeseg; %%
% nintn = par.geo.nintn;
% nnodes = par.geo.nnode;
% tns = (nis * nintn) + (nns * nnodes);
% tnf = tns - 1;
% 
% 
% % Indices along axon of nodes and internodes.
% % Indices of connections node-to-internode, internode-to-node and
% % internode-to-internode.
% nodes = reshape(bsxfun(@plus, (1:nns)', (0:nnodes-1)*(nns+nis)), 1, []);
% intn = 1 : tns;
% intn(nodes) = [];
% 
% nseg = nis + nns;
% seg = 0:(nintn-1);
% jpns = [nseg*seg+4; nseg*seg+5; seg*nseg+(nseg-3); seg*nseg+(nseg-2)];
% jpns = jpns(:)';
% 
% node2node = reshape(bsxfun(@plus, (1:(nns-1))', (0:nnodes-1)*(nns+nis)), 1, []);
% node2intn = nns:(nns+nis):tnf;
% intn2node = (nns+nis):(nns+nis):tnf;
% intn2intn = 1 : tnf;
% intn2intn([node2node, node2intn, intn2node]) = [];

% %2jpn jpn2 juxtaparanode
% intn2jpn = sort([3:(nns+nis):tnf,(nns+nis-4):(nns+nis):tnf]);
% jpn2intn = sort([5:(nns+nis):tnf,(nns+nis-2):(nns+nis):tnf]);


%% 
% figure
% hold on
% xline(jpns,'r', 'HandleVisibility','off')
% plot(1:tns+2, V2(:,2),'k', 'DisplayName', 'voltage')
% yline(Vlowerlimit,'--k','DisplayName','VLowLim')
% yline(Vupperlimit,'--k','DisplayName','VUpLim')
% xlabel('segment')
% ylabel('V2')
% legend
% hold off

