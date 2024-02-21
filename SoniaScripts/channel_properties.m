clear all
close all
clc

par = Cullen2018CortexAxonJPNlocalized_MTR;
nd= 15;

par.sim.dt.value = 1;
par.sim.tmax.value = 10;
[MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;


%% FIGURES (faster without the simulation above)
figure(1)
subplot(211)
hold on
plot(TIME_VECTOR,MEMBRANE_POTENTIAL(:,nd), 'k', 'LineWidth',1.5);
title('voltage')
xlabel('time (ms)')
ylabel('voltage')
hold off

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
        figure(1)
        subplot(212)
        hold on
        plot(TIME_VECTOR,gates{i}{j},'DisplayName', par.channels(i).channame, 'LineWidth', 1.5)
    end
end
hold off
legend
title('gating var')
xlabel('time (ms)')

%% K+ gates specifically
V = linspace(-100,100,length(TIME_VECTOR));

i = 4; 
g4 = nan(length(V),2);
a = nan(length(V), 1);
b = nan(length(V), 1);


for j=1:par.channels(i).gates.number
    for v=1:length(V)
        a(v) = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
        b(v) = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
        g4(v,j)= a(v)/(a(v)+b(v));
    end
    figure(2)
    subplot(321)
    hold on
    plot(V,a,'LineWidth',1.5, 'DisplayName',par.channels(i).gates.label{j})
    hold off
    legend 
    xlabel('mV')
    title('alpha')
    
    if j<2
        figure(2)
        subplot(323)
        plot(V,b,'LineWidth',1.5)
        xlabel('mV')
        title('beta - l')

        a1=a;
    else
        figure(2)
        subplot(325)
        plot(V,b, 'r','LineWidth',1.5)
        xlabel('mV')
        title('beta - r')
    end

    figure(2)
    subplot(322)
    hold on
    plot(V, g4(:,j),'LineWidth', 1.5, 'DisplayName', par.channels(i).gates.label{j})
    hold off
    legend
    xlabel('mV')
    title('a/(a+b)')
end
%
figure(2)
subplot(324)
plot(V,a1.*a, 'LineWidth', 1.5)
xlabel('mV')
title('alpha(l)*alpha(r)')

