function par = parameters_Richardson



%% simulation

% VALUES
par.sim.dt.value =                                                      1;
par.sim.tmax.value =                                                    1;
par.sim.temp =                                                          37;

par.stim.amp.value =                                                    192.9150825356307;
par.stim.dur.value =                                                    1;
par.stim.location =                                                     1;


par.geo.nnode =                                                         21;
par.geo.nintn =                                                         par.geo.nnode - 1;
par.geo.nnodeseg =                                                      1;
par.geo.nintseg =                                                       10;
par.geo.totalNumberSegments =                                           par.geo.nnode*par.geo.nnodeseg + par.geo.nintn*par.geo.nintseg;


% UNITS
par.sim.dt.units =                                                      {1,'us',1};
par.sim.tmax.units =                                                    {1,'ms',1};

par.stim.amp.units =                                                    {2,'mA', 'cm', [1, -2]};
par.stim.dur.units =                                                    {1,'us',1};


%% information about indexing the axon


par.geo.nodeSegments = mat2cell( ...
                                bsxfun( ...
                                        @plus, ...
                                        (1:par.geo.nnodeseg), ...
                                        (par.geo.nnodeseg + par.geo.nintseg)*(0:par.geo.nintn)' ...
                                      ), ...
                                ones(1, par.geo.nnode), ...
                                par.geo.nnodeseg);

par.geo.internodeSegments = mat2cell( ...
                                     bsxfun( ...
                                            @plus, ...
                                            (1:par.geo.nintseg), ...
                                            par.geo.nnodeseg + (par.geo.nnodeseg + par.geo.nintseg)*(0:par.geo.nintn-1)'), ...
                                     ones(1, par.geo.nintn), ...
                                     par.geo.nintseg);



%% geometry

% VALUES
par.node.geo.diam.value.ref =                                           3.3;
par.node.geo.length.value.ref =                                         1;

par.intn.geo.diam.value.ref=                                            6;
par.intn.geo.length.value.ref=                                          1149;

par.myel.geo.peri.value.ref =                                           10;
par.myel.geo.period.value =                                             16.652719665272;
par.myel.geo.numlamellae.value.ref =                                    120;


% UNITS
par.node.geo.diam.value.vec =                                           par.node.geo.diam.value.ref * ones(par.geo.nnode,1);
par.node.geo.diam.units =                                               {1,'um',1};
par.node.geo.length.value.vec =                                         par.node.geo.length.value.ref * ones(par.geo.nnode,1);
par.node.geo.length.units =                                             {1,'um',1};


par.intn.geo.diam.value.vec =                                           par.intn.geo.diam.value.ref * ones(par.geo.nintn,1);
par.intn.geo.diam.units =                                               {1,'um',1};
par.intn.geo.length.value.vec =                                         par.intn.geo.length.value.ref * ones(par.geo.nintn, 1);
par.intn.geo.length.units =                                             {1,'um',1};


par.node.seg.geo.length.value.vec =                                     repmat(par.node.geo.length.value.vec ./ par.geo.nnodeseg, 1, par.geo.nnodeseg);
par.node.seg.geo.length.units =                                         {1,'um',1};
par.node.seg.geo.diam.value.vec =                                       repmat(par.node.geo.diam.value.vec, 1, par.geo.nnodeseg);
par.node.seg.geo.diam.units =                                           {1,'um',1};


par.intn.seg.geo.length.value.vec =                                     repmat(par.intn.geo.length.value.vec ./ par.geo.nintseg, 1, par.geo.nintseg);
par.intn.seg.geo.length.units =                                         {1,'um',1};
par.intn.seg.geo.diam.value.vec =                                       repmat(par.intn.geo.diam.value.vec, 1, par.geo.nintseg);
par.intn.seg.geo.diam.units =                                           {1,'um',1};


par.myel.geo.peri.value.vec =                                           par.myel.geo.peri.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.geo.peri.units =                                               {1,'nm',1};
par.myel.geo.numlamellae.value.vec =                                    par.myel.geo.numlamellae.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.geo.period.units =                                             {1,'nm',1};





%% passive electrical

% VALUES
par.elec.pas.vrest.value.ref =                                          -82;

par.node.elec.pas.cap.value.ref =                                       2;
par.node.elec.pas.cond.value.ref =                                      80;
par.node.elec.pas.axres.value.ref =                                     70;

par.intn.elec.pas.cap.value.ref =                                       1;
par.intn.elec.pas.cond.value.ref =                                      0.1;
par.intn.elec.pas.axres.value.ref =                                     70;

par.myel.elec.pas.cap.value.ref =                                       0.1;
par.myel.elec.pas.cond.value.ref =                                      1;
par.myel.elec.pas.axres.value.ref =                                     70;


% UNITS
par.elec.pas.vrest.value.vec =                                          par.elec.pas.vrest.value.ref * ones(par.geo.totalNumberSegments, 1);
par.elec.pas.vrest.units =                                              {1, 'mV', 1};

par.node.elec.pas.cap.value.vec =                                       par.node.elec.pas.cap.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.cap.units =                                           {2,'uF','cm',[1 -2]};
par.node.elec.pas.cond.value.vec =                                      par.node.elec.pas.cond.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.cond.units =                                          {2,'mS','cm',[1 -2]};

par.node.elec.pas.axres.value.vec =                                     par.node.elec.pas.axres.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.axres.units =                                         {2,' O','cm',[1 1]};

par.intn.elec.pas.cap.value.vec =                                       par.intn.elec.pas.cap.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.cap.units =                                           {2,'uF','cm',[1 -2]};
par.intn.elec.pas.cond.value.vec =                                      par.intn.elec.pas.cond.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.cond.units =                                          {2,'mS','cm',[1 -2]};
par.intn.elec.pas.axres.value.vec =                                     par.intn.elec.pas.axres.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.axres.units =                                         {2,' O','cm',[1 1]};

par.myel.elec.pas.cap.value.vec =                                       par.myel.elec.pas.cap.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.cap.units =                                           {2,'uF','cm',[1,-2]};
par.myel.elec.pas.cond.value.vec =                                      par.myel.elec.pas.cond.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.cond.units =                                          {2,'mS','cm',[1,-2]};
par.myel.elec.pas.axres.value.vec =                                     par.myel.elec.pas.axres.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.axres.units =                                         {2,' O','cm',[1 1]};




%% active electrical

par.channels(1) =                                   RichardsonMcIntyreGrill2000FastNa;
par.channels(2) =                                   RichardsonMcIntyreGrill2000PersistentNa;
par.channels(3) =                                   RichardsonMcIntyreGrill2000SlowK;
par.channels(4) =                                   RichardsonMcIntyreGrill2000SlowK;
par.channels(4).cond.value = 0 * par.channels(4).cond.value;

par.channels(1).location =                          [par.geo.nodeSegments{:}]';
par.channels(2).location =                          [par.geo.nodeSegments{:}]';
par.channels(3).location =                          [par.geo.nodeSegments{:}]';
par.channels(4).location =                          [par.geo.internodeSegments{:}]';


for i = 1 : length(par.channels)
    par.channels(i).cond.value =                    par.channels(i).cond.value * ones(size(par.channels(i).location));
end

if isempty(par.channels)
    par.active.segments =                           [];
else
    par.active.segments =                           unique(cat(1, par.channels(:).location));
end


for i = 1 : length(par.channels)
    par.channels(i).locationInActive = [];
    for j = 1 : length(par.channels(i).location)
        assert(any(par.active.segments == par.channels(i).location(j)))
        assert(sum(par.active.segments == par.channels(i).location(j)) < 2)
        par.channels(i).locationInActive(j) =       find(par.active.segments == par.channels(i).location(j));
    end
end


gLN = unitsabs(par.node.elec.pas.cond.units) * par.node.elec.pas.cond.value.vec;
gLI = unitsabs(par.intn.elec.pas.cond.units) * par.intn.elec.pas.cond.value.vec;
par.elec.pas.cond.value.vec = [reshape([gLN(1:end-1, :), gLI]', [], 1); gLN(end, :)'];
par.elec.pas.cond.units = {2, 'mS', 'mm', [1, -2]};

clear gLN gLI

par.elec.pas.erev.value.vec =                                       updateleakerev(par);
par.elec.pas.erev.units =                                           {1,'mV',1};

