function par = Cullen2018CortexAxonJPNlocalized_MTR(jpn_cond)

% Initialize all parameters.
par =                                                                   GenerateEmptyParameterStructure();

% Simulation parameters
par.sim.temp =                                                          21;
par.sim.dt.value =                                                      0.1;
par.sim.dt.units =                                                      {1, 'us', 1};
par.sim.tmax.value =                                                    5;
par.sim.tmax.units =                                                    {1, 'ms', 1};

% Current stimulation
% Stimulation amplitude.
par.stim.amp.value =                                                    0.5;
par.stim.amp.units =                                                    {1, 'nA', 'cm', [1, -2]};
% Stimulation duration.
par.stim.dur.value =                                                    10;
par.stim.dur.units =                                                    {1, 'us', 1};
par.stim.location =                                                     1;

% Number of nodes.
par.geo.nnode =                                                         51;
% Number of internodes.
par.geo.nintn =                                                         par.geo.nnode - 1;
% Number of segments per node.
par.geo.nnodeseg =                                                      1;
% Number of segments per internode.
par.geo.nintseg =                                                       52;

% Number of juxtaparanodes.
par.geo.njpn =                                                          par.geo.nintn*2;
% Number of segments per juxtaparanode region
par.geo.jpnspacing =                                                    2;
par.geo.njpnseg =                                                       2;

% Total number of segments
par.geo.totalNumberSegments =                                           par.geo.nnode*par.geo.nnodeseg + par.geo.nintn*par.geo.nintseg;

% Information about indexing of the axon
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

nseg = par.geo.nintseg + par.geo.nnodeseg;
seg = 0:(par.geo.nintn-1);

jpns = [];
for i = 0:par.geo.njpnseg-1
    jpns = [jpns; nseg*seg + 2 + par.geo.jpnspacing + i;
        seg*nseg + (nseg - par.geo.jpnspacing - i)];
end
jpns = jpns';
% jpns = [nseg*seg+4; nseg*seg+5; seg*nseg+(nseg-3); seg*nseg+(nseg-2)]';
par.geo.juxtaparanodeSegments = mat2cell(jpns, ones(1, par.geo.nintn), par.geo.njpnseg*2);

% Node geometry
% Node diameter.
par.node.geo.diam.value.ref =                                           0.5894;
par.node.geo.diam.value.vec =                                           par.node.geo.diam.value.ref * ones(par.geo.nnode, 1);
par.node.geo.diam.units =                                               {1, 'um', 1};

% Node length.
par.node.geo.length.value.ref =                                         0.8364;
par.node.geo.length.value.vec =                                         par.node.geo.length.value.ref * ones(par.geo.nnode, 1);
par.node.geo.length.units =                                             {1, 'um', 1};

% Node segment diameter.
par.node.seg.geo.diam.value.ref =                                       par.node.geo.diam.value.ref;
par.node.seg.geo.diam.value.vec =                                       repmat(par.node.geo.diam.value.vec, 1, par.geo.nnodeseg);
par.node.seg.geo.diam.units =                                           {1, 'um', 1};

% Node segment length.
par.node.seg.geo.length.value.ref =                                     par.node.geo.length.value.ref;
par.node.seg.geo.length.value.vec =                                     repmat(par.node.geo.length.value.vec / par.geo.nnodeseg, 1, par.geo.nnodeseg);
par.node.seg.geo.length.units =                                         {1, 'um', 1};

% Internode geometry
% Internode axon diameter.
par.intn.geo.diam.value.ref =                                           0.5894;
par.intn.geo.diam.value.vec =                                           par.intn.geo.diam.value.ref * ones(par.geo.nintn, 1);
par.intn.geo.diam.units =                                               {1, 'um', 1};
% Internode length.
par.intn.geo.length.value.ref =                                         50.32;%48.32;
par.intn.geo.length.value.vec =                                         par.intn.geo.length.value.ref * ones(par.geo.nintn, 1);
par.intn.geo.length.units=                                              {1, 'um', 1};
% Internode segment length.
par.intn.seg.geo.length.value.ref =                                     par.intn.geo.length.value.ref / par.geo.nintseg;
par.intn.seg.geo.length.value.vec =                                     repmat(par.intn.geo.length.value.vec / par.geo.nintseg, 1, par.geo.nintseg);
par.intn.seg.geo.length.units =                                         {1, 'um', 1};
% Internode segment diameter (=internode diameter).
par.intn.seg.geo.diam.value.ref =                                       par.intn.geo.diam.value.ref;
par.intn.seg.geo.diam.value.vec =                                       repmat(par.intn.geo.diam.value.vec, 1, par.geo.nintseg);
par.intn.seg.geo.diam.units =                                           {1, 'um', 1};

% General electrical
% Resting membrane potential.
par.elec.pas.vrest.value.ref =                                          -72;
par.elec.pas.vrest.value.vec =                                          par.elec.pas.vrest.value.ref * ones(par.geo.totalNumberSegments, 1);
par.elec.pas.vrest.units =                                              {1, 'mV', 1};


% Node leak reversal potential.
% par.node.elec.pas.leak.erev.value.ref =                                 -84;
% par.node.elec.pas.leak.erev.value.vec =                                 par.node.elec.pas.leak.erev.value.ref * ones(par.geo.nnode, par.geo.nnodeseg);
% par.node.elec.pas.leak.erev.units =                                     {1, 'mV', 1};

% Node leak conductance - adjusted to set resting membrane potential.
par.node.elec.pas.cond.value.ref =                                      0.3211; % same as Richardson - from correspondence -- don't actually know what it should be
par.node.elec.pas.cond.value.vec =                                      par.node.elec.pas.cond.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.cond.units =                                          {2, 'mS', 'mm', [1, -2]};

% fileStr = {'McIntyre2002FastNa.mat', 'McIntyre2002PersistentNa.mat', 'McIntyre2002SlowK.mat'};
% for fileIdx = 1 : length(fileStr)
%     par =                                                               AddActiveChannel(par, fileStr{fileIdx});
% end
% par =                                                                   CalculateLeakConductance(par);

% Node axial resistivity.
par.node.elec.pas.axres.value.ref =                                     0.7;
par.node.elec.pas.axres.value.vec =                                     par.node.elec.pas.axres.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.axres.units =                                         {2, ' O', ' m', [1, 1]};

% Node membrane capacitance.
par.node.elec.pas.cap.value.ref =                                       0.9;
par.node.elec.pas.cap.value.vec =                                       par.node.elec.pas.cap.value.ref * ones(par.geo.nnode,par.geo.nnodeseg);
par.node.elec.pas.cap.units =                                           {2, 'uF', 'cm', [1, -2]};

% Myelin membrane capacitance.
par.myel.elec.pas.cap.value.ref =                                       0.9;
par.myel.elec.pas.cap.value.vec =                                       par.myel.elec.pas.cap.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.cap.units =                                           {2, 'uF', 'cm', [1, -2]};

% Myelin membrane conductance.
par.myel.elec.pas.cond.value.ref =                                      1;
par.myel.elec.pas.cond.value.vec =                                      par.myel.elec.pas.cond.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.cond.units =                                          {2, 'mS', 'cm', [1, -2]};

% Internode axon membrane capacitance.
par.intn.elec.pas.cap.value.ref =                                       0.9;
par.intn.elec.pas.cap.value.vec =                                       par.intn.elec.pas.cap.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.cap.units =                                           {2, 'uF', 'cm', [1, -2]};

% Internode axon membrane conductance.
par.intn.elec.pas.cond.value.ref =                                      0.1;
par.intn.elec.pas.cond.value.vec =                                      par.intn.elec.pas.cond.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.cond.units =                                          {2, 'mS', 'cm', [1, -2]}; 

% Internode axon membrane resistivity.
par.intn.elec.pas.axres.value.ref =                                     0.7;
par.intn.elec.pas.axres.value.vec =                                     par.intn.elec.pas.axres.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.intn.elec.pas.axres.units =                                         {2, ' O', ' m', [1, 1]};

% Periaxonal space resistivity.
par.myel.elec.pas.axres.value.ref =                                     0.7;
par.myel.elec.pas.axres.value.vec =                                     par.myel.elec.pas.axres.value.ref * ones(par.geo.nintn, par.geo.nintseg);
par.myel.elec.pas.axres.units =                                         {2, ' O', ' m', [1, 1]};

% Periaxonal space width.
par.myel.geo.peri.value.ref =                                           6.477;
par.myel.geo.peri.value.vec =                                           par.myel.geo.peri.value.ref * ones(par.geo.nintn,par.geo.nintseg);
par.myel.geo.peri.units =                                               {1, 'nm', 1};

% Myelin wrap periodicity.
par.myel.geo.period.value =                                             16.2873;
par.myel.geo.period.units =                                             {1, 'nm', 1};

% g-ratio (internode axon diameter to internode outer diameter ratio)
par.myel.geo.gratio.value.ref =                                         0.724;
par.myel.geo.gratio.value.vec_ref =                                     par.myel.geo.gratio.value.ref * ones(par.geo.nintn, par.geo.nintseg);

% Set units of myelin width.
par.myel.geo.width.units =                                              {1, 'um', 1};

% Update number of myelin lamellae in separate function
par =                                                                   CalculateNumberOfMyelinLamellae(par, 'max');

% Restrict periaxonal space around the paranodes.
% We need to use 'min' to update the number of myelin lamellae, as reducing
% the periaxonal space increases the number of layers in the paranodes.
% See CALCULATENUMBEROFMYELINLAMELLAE.m for details.
% par =                                                                   UpdateInternodePeriaxonalSpaceWidth(par, par.myel.geo.peri.value.ref/2, [], [1, 2, 51, 52], 'min');


% active electrical - introduced from correspondence

par.channels(1) =                                   McIntyre2002FastNa_JPN;
par.channels(2) =                                   McIntyre2002PersistentNa_JPN;
par.channels(3) =                                   McIntyre2002SlowK_JPN;
par.channels(4) =                                   Kv11_JPN;
par.channels(5) =                                   Kv12_JPN;

if nargin > 0
    par.channels(4).cond.value = jpn_cond;
    par.channels(5).cond.value = jpn_cond;
end

par.channels(1).location =                          [par.geo.nodeSegments{:}]';
par.channels(2).location =                          [par.geo.nodeSegments{:}]';
par.channels(3).location =                          [par.geo.nodeSegments{:}]';
par.channels(4).location =                          [par.geo.juxtaparanodeSegments{:}]';
par.channels(5).location =                          [par.geo.juxtaparanodeSegments{:}]';


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


