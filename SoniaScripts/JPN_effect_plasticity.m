clear
clc
close all
% Initiate directory for saving data.
thisDirectory   = fileparts(mfilename('fullpath'));
saveDirectory   = fullfile(thisDirectory,'MTR_JPN_R4.Kv12_2.0S');
if ~isdir(saveDirectory)
    mkdir(saveDirectory)
end


%
cond_vals = [0 0.02 3]; % S/cm2
velocity = nan(length(cond_vals),2,4);

for g = 1:length(cond_vals)
    clear par1  par2  par_opt

    par1 = Cullen2018CortexAxonJPNlocalized_MTR(cond_vals(g));
    par2 = Cullen2018CortexAxonJPNlocalized_Kv12(cond_vals(g));
    par_opt = [par1; par2];
    for k = 1:2

        % Produce parameters for default cortex model.
        clear par;
        par = par_opt(k);
        par.sim.dt.value = 1;
        % Set temperature.
        par.sim.temp = 37;
        % Adjust simulation time.
        par.sim.tmax.value = 15;

        %% Run all simulations varying periaxonal space width.
        j = 1;
        for psw = [0:0.2:1.6 2:6 6.477 7:8 8.487 9 10:2:14 15 20]
            par.myel.geo.peri.value.ref = psw;
            par.myel.geo.peri.value.vec = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
            par.myel.geo.period.value   = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
            par                         = CalculateNumberOfMyelinLamellae(par, 'max');
            %         par                         = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
            [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par, fullfile(saveDirectory, ['MTR2024_psw_' num2str(psw) '_' num2str(par.sim.temp) 'C.mat']));
            velocity_psw(g,j,1)           = psw;
            if k == 1
                velocity_psw(g,j,2)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
            else
                velocity_psw(g,j,3)       = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);
            end
            j                           = j+1;
        end
        refresh;

        % reset model and temperature
        clear par;
        par = par_opt(k);
        par.sim.dt.value = 1;
        % Set temperature.
        par.sim.temp = 37;

        %% Run sham simulations.
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par, fullfile(saveDirectory, ['MTR2024_sham_' num2str(par.sim.temp) 'C.mat']));

        velocity(g,k,1) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

        refresh;

        % Reset model and temperature.
        clear par;
        par             = par_opt(k);
        par.sim.temp    = 37;
        par.sim.dt.value = 1;
        %% Run short node simulations.
        par.node.geo.length.value.ref       = 0.7735;
        par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
        par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
        par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
        par =                                 CalculateLeakConductanceJPN_MTR(par);
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par, fullfile(saveDirectory, ['MTR2024_shortnode_' num2str(par.sim.temp) 'C.mat']));

        velocity(g,k,2) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

        refresh;

        % Reset model and temperature.
        clear par;
        par             = par_opt(k);
        par.sim.temp    = 37;
        par.sim.dt.value = 1;
        %% Run alt. myelin simulations.
        par.myel.geo.gratio.value.ref       = 0.6888;
        par.myel.geo.gratio.value.vec_ref   = par.myel.geo.gratio.value.ref * ones(par.geo.nintn, par.geo.nintseg);
        par.myel.geo.peri.value.ref         = 8.487;
        par.myel.geo.peri.value.vec         = par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
        par.myel.geo.period.value           = 1000*(par.node.geo.diam.value.ref/par.myel.geo.gratio.value.ref-par.node.geo.diam.value.ref-2*par.myel.geo.peri.value.ref/1000)/(2*6.5);
        par                                 = CalculateNumberOfMyelinLamellae(par, 'max');
        %     par                                 = UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par, fullfile(saveDirectory, ['MTR2024_altmyelin_' num2str(par.sim.temp) 'C.mat']));

        velocity(g,k,3) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

        %% Run iTBS simulations.
        par.node.geo.length.value.ref       = 0.7735;
        par.node.geo.length.value.vec       = par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
        par.node.seg.geo.length.value.ref   = par.node.geo.length.value.ref;
        par.node.seg.geo.length.value.vec   = repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
        par =                                 CalculateLeakConductanceJPN_MTR(par);
        [MEMBRANE_POTENTIAL, INTERNODE_LENGTH, TIME_VECTOR] = ModelJPN_MTR(par, fullfile(saveDirectory, ['MTR2024_iTBS_' num2str(par.sim.temp) 'C.mat']));

        velocity(g,k,4) = velocities(MEMBRANE_POTENTIAL, INTERNODE_LENGTH, par.sim.dt.value*simunits(par.sim.dt.units), [20 40]);

    end
end
save("velocity_JPN_effect.mat", "velocity")
save("velocity_psw_JPN_effect.mat", "velocity_psw")



%% Plotting
% velocity = load("velocity_JPN_effect.mat");
% velocity_psw = load("velocity__psw_JPN_effect.mat")
% cond_vals = [0 0.02 3]; % S/cm2
dif = nan(3,3,2);

for g = 2:length(cond_vals)
    velocity_0 = reshape(velocity(1,1,:), [1, 4]);
    velocity_g = [velocity_0; reshape(velocity(g,:,:), [2, 4])]';

    for i = 1:3
        dif(i,:,g-1) = 100*(velocity_g(i+1,:)./velocity_g(1,:)-1);
    end

    if g == 2
        subplot(2,2,1);
        hold on
        plot(velocity_psw(1,:,1),velocity_psw(1,:,3),'color',[0 0.4470 0.7410], LineWidth=1.2),
        plot(velocity_psw(g,:,1),velocity_psw(g,:,2),'-k', LineWidth=1.2), 
        plot(velocity_psw(g,:,1),velocity_psw(g,:,3),'-r', LineWidth=1.2);

        % xlabel('psw (nm)'),
        ylabel('CV(m/s)')
        legend('no JPN','Kv1.1', 'Kv1.2')
        title(['Conductance ' num2str(cond_vals(g)) ' S/cm^2'])
        xline(6.477,'HandleVisibility','off', LineWidth=1.1)

        subplot(2,2,2);
        hold on;
        b = bar(velocity_g,'FaceColor','flat');
        b(1).CData = [0 0.4470 0.7410];  
        b(2).CData = [0 0 0];            
        b(3).CData = [1 0 0];   
        for i = 2:4
            plot(i,1.95,'v','MarkerEdgeColor','k','MarkerFaceColor','k'); %1.17
            text(i,2.1,[num2str(100*(velocity_g(i,1)./velocity_g(1,1)-1),'%1.1f') '%'],'Color','k','HorizontalAlignment','center','FontSize',8); %1.185
        end
        xticks(1:1:4.5), xticklabels({'Control','Short node','Alt. myelin','iTBS'}),
        yticks(0:0.5:2.5)%, yticklabels({'','0.5''1','1.5'}),
        ylabel('CV (m/s)');
        ylim([0 2.3])

    else
        subplot(2,2,3);
        hold on
        plot(velocity_psw(1,:,1),velocity_psw(1,:,3),'color',[0 0.4470 0.7410], LineWidth=1.2),
        plot(velocity_psw(g,:,1),velocity_psw(g,:,2),'-k', LineWidth=1.2), 
        plot(velocity_psw(g,:,1),velocity_psw(g,:,3),'-r', LineWidth=1.2);

        xlabel('psw (nm)'), ylabel('CV(m/s)')
        legend('no JPN','Kv1.1', 'Kv1.2')
        title(['Conductance ' num2str(cond_vals(g)) ' S/cm^2'])
        xline(6.477,'HandleVisibility','off', LineWidth=1.1)

        subplot(2,2,4);
        hold on
        b = bar(velocity_g,'FaceColor','flat');
        b(1).CData = [0 0.4470 0.7410];  
        b(2).CData = [0 0 0];            
        b(3).CData = [1 0 0];   
        for i = 2:4
            plot(i,1.95,'v','MarkerEdgeColor','k','MarkerFaceColor','k'); %1.17
            % text(i,2.1,[num2str(100*(velocity_g(i,1)./velocity_g(1,1)-1),'%1.1f') '%'],'Color','k','HorizontalAlignment','center','FontSize',8); %1.185
        end
        xticks(1:1:4.5), xticklabels({'Control','Short node','Alt. myelin','iTBS'}),
        yticks(0:0.5:2.5)%, yticklabels({'','0.5''1','1.5'}),
        ylabel('CV (m/s)');
        ylim([0 2.3])
    end
end
