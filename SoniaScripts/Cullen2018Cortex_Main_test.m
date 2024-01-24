clear

% Initiate temperature.
temp            = 21;
k             = temp;

% Figure.
f               = figure;

% Produce parameters for default cortex model.
clear par;
par = Cullen2018CortexAxonJPN_MTR();
    
% Set temperature.
par.sim.temp = k;
    
% Adjust simulation time.
par.sim.tmax.value = 10;
    
%% Run controls for insets.
perisw = [0 20];
for i= 1:length(perisw)
    psw = perisw(i);
    par.myel.geo.peri.value.ref = psw;
    par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
    par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
    par                         = CalculateNumberOfMyelinLamellae(par, 'max');
    %         par                         = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
    [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par);

    subplot(1,2,i)
    plot(TIME_VECTOR,MEMBRANE_POTENTIAL,'-k');
    title(["psw: " num2str(psw)]);
end
% refresh;


