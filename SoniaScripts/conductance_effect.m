
%% ----------- AP SHAPE (large step sim) --------------

cond_vals = 0:0.05:2.5; %[0 0.05 0.5 1 1.75 2.5]; %S/cm2
psw = 6.477; %, 0:5:20]; %nm !!! changing this will mess up plotting
% spacing, n. jpn seg
config = [1 2];

% Color Scales
red = [1, 0, 0];
blue = [0, 0, 1];

colorMapLength_p = length(psw)-1;
color_scale_p = [linspace(red(1),blue(1),colorMapLength_p)', linspace(red(2),blue(2),colorMapLength_p)', linspace(red(3),blue(3),colorMapLength_p)'];
color_scale_p = [0, 0, 0; color_scale_p];

colorMapLength_g = length(cond_vals)-1;
color_scale_g = [linspace(red(1),blue(1),colorMapLength_g)', linspace(red(2),blue(2),colorMapLength_g)', linspace(red(3),blue(3),colorMapLength_g)'];
color_scale_g = [0, 0, 0; color_scale_g];


velocity = nan(length(cond_vals),length(psw),size(config,1));
tmax = 300;

% figure(1)
for c = 1:size(config,1)
    for p = 1:length(psw)
        for g = 1:length(cond_vals)
            clear par
            % parameters
            par = Cullen2018CortexAxonJPNlocalized_MTR(cond_vals(g));
                % geomery
            par.geo.jpnspacing = config(c,1);
            par.geo.njpnseg = config(c,2);
                % psw
            par.myel.geo.peri.value.ref = psw(p);
            par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
            par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
            par                         = CalculateNumberOfMyelinLamellae(par, 'max');
                % simulation parameters
            par.sim.temp        = 37;
            par.sim.dt.value    = 5;
            par.sim.tmax.value  = tmax; %800

            % simulation
            [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par);

            % plotting
            % figure(1)
            % % subplot(3, 2, p)
            % hold on
            % plot(TIME_VECTOR,MEMBRANE_POTENTIAL(:,3),'color', color_scale_g(g,:), 'DisplayName', num2str(cond_vals(g)));
            % hold off
            % title([num2str(psw(p)) ' nm'])
            % legend
            % xlabel('Time (ms)'), ylabel('Axon voltage (mV)');

            velocity(g,p,c) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
    end
    % check = 1;
end
%%
figure(2)
for c = 1:size(config,1)
    for p = 1:length(psw)
        % subplot(3, 2, p)
        hold on
        plot(cond_vals,velocity(:,p,c),'-k','DisplayName','Kv1.1');
        xline(0.05,'--r', 'DisplayName','hyperexcitability boundary')
        xline(1.6, '--r', DisplayName='hyperexcitability boundary')
        title([num2str(psw(p)) ' nm'])
        xlabel('conductance'), ylabel('velocity')
        ylim([1.4, 2.1])
        legend
        hold off

    end
end

%%  --------------- CV EFFECT ----------------------

cond_vals = 0:0.05:1.7; %S/cm2
psw = 6.477; %, 0:5:20]; %nm !!! changing this will mess up plotting

velocity = nan(length(cond_vals),length(psw),size(config,1));
tmax = 5;

% figure(1)
for c = 1:size(config,1)
    for p = 1:length(psw)
        for g = 1:length(cond_vals)
            clear par
            % parameters
            par = Cullen2018CortexAxonJPNlocalized_Kv12(cond_vals(g));
                % geomery
            par.geo.jpnspacing = config(c,1);
            par.geo.njpnseg = config(c,2);
                % psw
            par.myel.geo.peri.value.ref = psw(p);
            par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
            par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
            par                         = CalculateNumberOfMyelinLamellae(par, 'max');
                % simulation parameters
            par.sim.temp        = 37;
            par.sim.dt.value    = 1;
            par.sim.tmax.value  = tmax; %800

            % simulation
            [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par);
            velocity(g,p,c) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
        end
    end
    % check = 1;
end
figure(3)
for c = 1:size(config,1)
    for p = 1:length(psw)
        % subplot(3, 2, p)
        hold on
        plot(cond_vals,velocity(:,p,c),'-k');
        title([num2str(psw(p)) ' nm'])
        xlabel('conductance'), ylabel('velocity')
        legend
        hold off

    end
end
