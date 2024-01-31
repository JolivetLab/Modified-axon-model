function par = CalculateLeakConductanceJPN_MTR(par)
%CALCULATELEAKCONDUCTANCE - Calculate leak conductance of nodes.
%   par = CALCULATELEAKCONDUCTANCE(par)
%       Outputs:
%           par -           modified parameter structure
%       Inputs:
%           par -           parameter structure
%   
%   Calculates the leak conductance required to set the resting membrane
%   potential, given a collection of active conductances in the nodes.
%   
%   The `par' structure is updated within this function so all necessary
%   parameters must be set before it is called (including the desired leak
%   conductance units).

% Axon geometry
PERIAXONAL_SPACE = simunits(par.myel.geo.peri.units) * par.myel.geo.peri.value.vec;

% Node and internode radius.
RADIUS_NODE = simunits(par.node.geo.diam.units)*par.node.seg.geo.diam.value.vec/2;
RADIUS_INODE = simunits(par.intn.seg.geo.diam.units)*par.intn.seg.geo.diam.value.vec/2;

% Myelin radius (and number of lamellae/myelin wrap periodicity)
NUMBER_LAMELLAE = par.myel.geo.numlamellae.value.vec;
PERIODICITY = simunits(par.myel.geo.period.units)*par.myel.geo.period.value;
RADIUS_MYELIN = cell(1,par.geo.nintn);
for i=1:par.geo.nintn
    RADIUS_MYELIN{i} = repmat(RADIUS_INODE(i,:),2*NUMBER_LAMELLAE(i,1),1)...
                        +repmat(PERIAXONAL_SPACE(i,:),2*NUMBER_LAMELLAE(i,1),1)...
                        +repmat(((1:2*NUMBER_LAMELLAE(i))'-1)*(PERIODICITY/2),1,par.geo.nintseg);
end

% Node and internode length
LENGTH_NODE = simunits(par.node.geo.length.units)*par.node.seg.geo.length.value.vec;                        % nnodex1 vec
LENGTH_INODE = simunits(par.intn.seg.geo.length.units)*par.intn.seg.geo.length.value.vec;               % nintnxnintseg mat

% Node, internode and myelin membrane surface areas.
SURFACEAREA_NODE = 2*pi*RADIUS_NODE.*LENGTH_NODE;                                                       % nnodex1 vec
SURFACEAREA_INODE = 2*pi*RADIUS_INODE.*LENGTH_INODE;                                                    % nintnxnintseg mat
SURFACEAREA_MYELIN=cell(1,par.geo.nintn);
for i=1:par.geo.nintn
    SURFACEAREA_MYELIN{i} = 2*pi*RADIUS_MYELIN{i}.*repmat(LENGTH_INODE(i,:),2*NUMBER_LAMELLAE(i,1),1);  % nlamxnintseg mat
end

surfaceAreaAxolemma = [reshape([SURFACEAREA_NODE(1:end-1, :), SURFACEAREA_INODE]', [], 1); SURFACEAREA_NODE(end, :)'];
%%
nodes = reshape(bsxfun(@plus, (1:par.geo.nnodeseg)', (0:par.geo.nnode-1)*(par.geo.nnodeseg+par.geo.nintseg)), 1, []);
%%

% Resting membrane potential
vrest = simunits(par.elec.pas.vrest.units) * par.elec.pas.vrest.value.ref;

% Active conductances
actcond=cell(1,length(par.channels));
for i=1:length(par.channels)
    actcond{i}=simunits(par.channels(i).cond.units)*par.channels(i).cond.value.*surfaceAreaAxolemma(par.channels(i).location);
end

% Gating variable values at resting membrane potential.
gates = cell(1, length(par.channels));
for i = 1 : length(par.channels)
    gates{i} = cell(1, par.channels(i).gates.number);
    for j = 1 : par.channels(i).gates.number
        a = rateequation(vrest, par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.alpha.q10(j), par.channels(i).gates.alpha.equ{j});
        b = rateequation(vrest, par.sim.temp, par.channels(i).gates.temp, par.channels(i).gates.beta.q10(j), par.channels(i).gates.beta.equ{j});
        gates{i}{j}(1:length(par.channels(i).location), 1) = a/(a+b);
    end
end

% Sum of active channel currents at resting membrane potential.
activesum = zeros(size(par.active.segments)); %zeros(par.geo.nnode, par.geo.nnodeseg);
for j = 1 : length(par.channels)
    tempprod = ones(length(par.channels(j).location),1); %ones(par.geo.nnode, par.geo.nnodeseg); 
    for k = 1 : par.channels(j).gates.number %par.node.elec.act(j).gates.number
        tempprod = tempprod.*(gates{j}{k}.^par.channels(j).gates.numbereach(k));
    end
    activesum(par.channels(j).locationInActive) = activesum(par.channels(j).locationInActive) + actcond{j} .* tempprod * (vrest - simunits(par.channels(j).erev.units) * par.channels(j).erev.value);
end

% Leak conductances to balance the active channel currents.
par.node.elec.pas.cond.value.vec = ...
    -activesum(nodes) * fromsimunits(par.node.elec.pas.cond.units) ./ (vrest - simunits(par.elec.pas.erev.units)*par.elec.pas.erev.value.vec(nodes));
par.node.elec.pas.cond.value.ref = mode(par.node.elec.pas.cond.value.vec(:));