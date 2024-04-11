clear all
close all
clc

for c = 1:2

    if c == 1
        par = Cullen2018CortexAxonJPNlocalized_MTR(0);
        par.sim.temp = 24;
    else
        par = Cullen2018CortexAxonJPNlocalized_Kv12(0);
        par.sim.temp = 20;
    end

    nd= 15;

    par.sim.dt.value = 1;
    par.sim.tmax.value = 15;
    [MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

    dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;


    %% FIGURES (faster without the simulation above)
    % par = Cullen2018CortexAxonJPNlocalized_dimer(0);
    if c == 1
        figure(1)
        subplot(4,3,[1:2,4:5])
        hold on
        plot(TIME_VECTOR,MEMBRANE_POTENTIAL(:,nd), 'k', 'LineWidth',1.5);
        set(gca, 'XTick', []);
        ylabel('Axon Voltage (mV)')
        hold off

        V = MEMBRANE_POTENTIAL(:,nd);
        vrest = simunits(par.elec.pas.vrest.units)*par.elec.pas.vrest.value.ref;
        gates=cell(1,length(par.channels));

        for i=3:length(par.channels)
            gates{i}=cell(1,par.channels(i).gates.number);
            for j=1:par.channels(i).gates.number
                for v=1:length(V)
                    a = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
                    b = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
                    gates{i}{j}((v), 1) = a/(a+b);
                end
                figure(1)
                subplot(4,3,[7:8,10:11])
                hold on
                plot(TIME_VECTOR,gates{i}{j},'DisplayName', par.channels(i).channame, 'LineWidth', 1.5)
            end
        end
        hold off
        legend
        leg = legend('show'); % Just one output here
        title(leg,'Channel')
        xlabel('Time (ms)'), ylabel('Gating Variable')
    else
        i=length(par.channels);
        gates{i}=cell(1,par.channels(i).gates.number);
        for j=1:par.channels(i).gates.number
            for v=1:length(V)
                a = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
                b = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
                gates{i}{j}((v), 1) = a/(a+b);
            end
            figure(1)
            subplot(4,3,[7:8,10:11])
            hold on
            plot(TIME_VECTOR,gates{i}{j},'DisplayName', par.channels(i).channame, 'LineWidth', 1.5)
        end
    end

    % K+ gates specifically
    V = linspace(-100,100,length(TIME_VECTOR));
    
    i = length(par.channels);
    g4 = nan(length(V),2);
    inf = nan(length(V), 1);
    tau = nan(length(V), 1);

    for j=1:par.channels(i).gates.number
        for v=1:length(V)
            inf(v) = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.inf.equ{j});
            tau(v) = rateequation(V(v), par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.tau.equ{j});
            g4(v,j)= inf(v)/(inf(v)+tau(v));
        end

        if j<2
            figure(1)
            subplot(4,3,3)
            hold on
            plot(V,inf,'DisplayName', par.channels(i).channame,'LineWidth',1.5)
            hold off
            set(gca, 'XTick', []);
            leg = legend('show'); % Just one output here
            leg.Location = 'southeast';
            title(leg,'Channel')
            title('m_{inf}')

            figure(1)
            subplot(4,3,9)
            hold on
            plot(V,tau,'LineWidth',1.5)
            hold off
            set(gca, 'XTick', []);
            title('\tau_m')

            inf1=inf;
        else
            figure(1)
            subplot(4,3,6)
            hold on
            plot(V,inf,'LineWidth',1.5)
            hold off
            set(gca, 'XTick', []);
            title('h_{inf}')

            figure(1)
            subplot(4,3,12)
            hold on
            plot(V, tau,'LineWidth',1.5)
            hold off
            xlabel('mV')
            title('\tau_h')
        end
    end
end
%
% figure(2)
% % subplot(325)
% plot(V,inf1.*inf, 'LineWidth', 1.5)
% xlabel('mV')
% title('inf(m)*inf(h)')

