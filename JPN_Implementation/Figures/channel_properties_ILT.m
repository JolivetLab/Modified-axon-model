clear all
close all
clc

clr = ["black", "red", "#77AC30"];

for c = 1:3

    if c == 1
        par = Cullen2018CortexAxonJPNlocalized_MTR(0);
        % par.sim.temp = 24;
        nd = 3;

        par.sim.temp = 24;
        par.sim.dt.value = 1;
        par.sim.tmax.value = 15;
        [MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

        dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
    elseif c == 2
        par = Cullen2018CortexAxonJPNlocalized_Kv12(0);
        % par.sim.temp = 20;
        nd = 3;

        par.sim.temp = 20;
        par.sim.dt.value = 1;
        par.sim.tmax.value = 15;
        [MEMBRANE_POTENTIAL, ~ , TIME_VECTOR] = ModelJPN_MTR(par);

        dt = unitsabs(par.sim.dt.units) * par.sim.dt.value;
        
    else
        par = Cullen2018CortexAxonJPNlocalized_KILT(0);
        par.sim.temp = 22;
    end


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
            subplot(321)
            hold on
            plot(V,inf,'Color', clr(c), 'DisplayName', par.channels(i).channame,'LineWidth',1.1)
            hold off
            % set(gca, 'XTick', []);
            % leg = legend('show'); % Just one output here
            % leg.Location = 'southeast';
            % title(leg,'Channel')
            title('m_{inf}')
            inf1=inf;

            figure(1)
            subplot(323)
            hold on
            plot(V,tau,'Color', clr(c), 'DisplayName', par.channels(i).channame,'LineWidth',1.1)
            hold off
            % set(gca, 'XTick', []);
            title('\tau_m (ms)')
            
        else
            figure(1)
            subplot(322)
            hold on
            plot(V,inf,'Color', clr(c), 'DisplayName', par.channels(i).channame,'LineWidth',1.1)
            hold off
            % set(gca, 'XTick', []);
            title('h_{inf}')   
            
            figure(1)
            subplot(3,2,5:6)
            hold on
            plot(V,inf1.*inf,'Color', clr(c), 'DisplayName', par.channels(i).channame,'LineWidth',1.1)
            xlabel('mV')
            title('m_{inf}*h_{inf}')
            hold off
            leg = legend('show'); % Just one output here
            % leg.Location = 'east';
            title(leg,'Channel')

            figure(1)
            subplot(324)
            hold on
            plot(V, tau,'Color', clr(c), 'DisplayName', par.channels(i).channame,'LineWidth',1.1)
            hold off
            xlabel('mV')
            title('\tau_h (ms)')
            % 
        end

    end

end


