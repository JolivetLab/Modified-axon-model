clear all
close all
clc

par = Cullen2018CortexAxonJPNlocalized_MTR;
nd= 15;

par.sim.dt.value = 1;
par.sim.tmax.value = 10;
[MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;

figure
subplot(211)
hold on
plot(TIME_VECTOR,MEMBRANE_POTENTIAL(:,nd));
title('voltage')
xlabel('time (ms)')
ylabel('voltage')
hold off
%%
V = MEMBRANE_POTENTIAL(:,nd);

vrest = simunits(par.elec.pas.vrest.units)*par.elec.pas.vrest.value.ref;

gates=cell(1,length(par.channels));

for i=1:length(par.channels)
    gates{i}=cell(1,par.channels(i).gates.number);
    for j=1:par.channels(i).gates.number
        for v=1:length(V)
            a = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
            b = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
            gates{i}{j}((v), 1) = a/(a+b);
        end
    end
end

subplot(212)
hold on
for i=1:length(par.channels)
    for j = 1:par.channels(i).gates.number
        plot(TIME_VECTOR,gates{i}{j},'DisplayName', par.channels(i).channame)
    end
end
hold off
legend
title('gating var')
xlabel('time (ms)')

%% K+ gates specifically
V = linspace(-120,50,length(TIME_VECTOR));

i = 4; 
j = 1; 
g4 = nan(length(V),1);
for v=1:length(V)
    a = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
    b = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
    g4(v) = a/(a+b);
end

figure
plot(V,g4)
xlabel('voltage')
ylabel('gating var')

