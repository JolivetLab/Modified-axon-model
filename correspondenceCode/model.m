function Vsave = model(par, varargin)


if nargin == 1
    options.segmentsToSave = [par.geo.nodeSegments{:}]';
    options.saveWithinModel = false;
    options.measureEnergy = false;
else
    options = varargin{1};
end

if options.saveWithinModel
    fid = fopen([options.fileName, '.bin'], 'w');
%     fidINa = fopen([options.fileName, 'Na.bin'], 'w');
%     fidIL = fopen([options.fileName, 'Leak.bin'], 'w');
end





fprintf('loading parameters...\n');
tic;


nis = par.geo.nintseg;
nns = par.geo.nnodeseg;
nintn = par.geo.nintn;
nnodes = par.geo.nnode;
tns = par.geo.totalNumberSegments;
tnf = tns - 1;
nnis = nns+nis;

nodes = [par.geo.nodeSegments{:}];
intn = [par.geo.internodeSegments{:}];

node2node = reshape((repmat(1:nns-1, nnodes, 1) + repmat(nnis*(0:nnodes-1)', 1, nns-1))', [], 1)';
node2intn = nns : nnis : tnf;
intn2node = nnis : nnis : tnf;
intn2intn = 1 : tnf;
intn2intn([node2node, node2intn, intn2node])=[];



dt = unitsabs(par.sim.dt.units)*par.sim.dt.value;
tmax = unitsabs(par.sim.tmax.units)*par.sim.tmax.value;
T = round(tmax/dt);





% geometry
PERIAXONAL_SPACE = unitsabs(par.myel.geo.peri.units) * par.myel.geo.peri.value.vec;
NUMBER_LAMELLAE = par.myel.geo.numlamellae.value.vec;
PERIODICITY = unitsabs(par.myel.geo.period.units)*par.myel.geo.period.value;


RADIUS_NODE = unitsabs(par.node.seg.geo.diam.units)*par.node.seg.geo.diam.value.vec/2;
RADIUS_INODE = unitsabs(par.intn.seg.geo.diam.units)*par.intn.seg.geo.diam.value.vec/2;
RADIUS_MYELIN = cell(1, nintn);
for i = 1 : nintn
    RADIUS_MYELIN{i}=repmat(RADIUS_INODE(i,:),2*NUMBER_LAMELLAE(i,1),1)...
                        +repmat(PERIAXONAL_SPACE(i,:),2*NUMBER_LAMELLAE(i,1),1)...
                        +repmat(((1:2*NUMBER_LAMELLAE(i))'-1)*(PERIODICITY/2),1,nis);
end

LENGTH_NODE = unitsabs(par.node.seg.geo.length.units)*par.node.seg.geo.length.value.vec;                            % nnodex1 vec
LENGTH_INODE = unitsabs(par.intn.seg.geo.length.units)*par.intn.seg.geo.length.value.vec;                   % nintnxnintseg mat


SURFACEAREA_NODE = 2*pi*RADIUS_NODE.*LENGTH_NODE;                                                    % nnode x nnodeseg vec
SURFACEAREA_INODE = 2*pi*RADIUS_INODE.*LENGTH_INODE;                                                 % nintn x nintseg mat
SURFACEAREA_MYELIN=cell(1,nintn);
for i=1:nintn
    SURFACEAREA_MYELIN{i} = 2*pi*RADIUS_MYELIN{i}.*repmat(LENGTH_INODE(i,:),2*NUMBER_LAMELLAE(i,1),1); % nlamxnintseg mat
end

surfaceAreaAxolemma = [reshape([SURFACEAREA_NODE(1:end-1, :), SURFACEAREA_INODE]', [], 1); SURFACEAREA_NODE(end, :)'];

AREA_NODE = pi*RADIUS_NODE.^2;          % nnodexnnodeseg vec
AREA_INODE = pi*RADIUS_INODE.^2;        % nintnxnintseg mat
AREA_PERI = pi*((RADIUS_INODE+PERIAXONAL_SPACE).^2-RADIUS_INODE.^2); % nintnxnintseg mat



Is = unitsabs(par.stim.amp.units) * par.stim.amp.value;
dur = unitsabs(par.stim.dur.units) * par.stim.dur.value;
DUR = ceil(dur/dt);

Istim = zeros(tns, 2, 2);
Istim([par.stim.location], 1, 1) = ...
    surfaceAreaAxolemma([par.stim.location])*Is;




% resting membrane potential
vrest = unitsabs(par.elec.pas.vrest.units)*par.elec.pas.vrest.value.vec;


% capacitance for each segment
C_NODE = unitsabs(par.node.elec.pas.cap.units)*par.node.elec.pas.cap.value.vec.*SURFACEAREA_NODE;    % nnode x nnodeseg vec
Cmy_INODE=cell(1,nintn);
CMY_INODE=zeros(nintn,nis);
for i = 1 : nintn
    Cmy_INODE{i} = unitsabs(par.myel.elec.pas.cap.units)*bsxfun(@times, par.myel.elec.pas.cap.value.vec(i, :), SURFACEAREA_MYELIN{i}); % nlamxnintseg mat
    CMY_INODE(i,:) = 1./(sum(1./Cmy_INODE{i},1));                                                         % nintnxnintseg mat
end   
CI_INODE = unitsabs(par.intn.elec.pas.cap.units)*par.intn.elec.pas.cap.value.vec.*SURFACEAREA_INODE;         % nintnxnintseg mat

cap=zeros(tns,2);
cap(nodes,1)=reshape(C_NODE', [], 1);
cap(intn,:)=[reshape(CI_INODE', nis*nintn, 1), reshape(CMY_INODE', nis*nintn, 1)];
clear C_NODE Cmy_INODE CMY_INODE CI_INODE





% leak conductance for each segment
LEAK_NODE = unitsabs(par.node.elec.pas.cond.units)*par.node.elec.pas.cond.value.vec.*SURFACEAREA_NODE;
LEAKmy_INODE=cell(1,nintn);
LEAKMY_INODE=zeros(nintn,nis);
for i = 1 : nintn
    LEAKmy_INODE{i} = unitsabs(par.myel.elec.pas.cond.units)*bsxfun(@times, par.myel.elec.pas.cond.value.vec(i, :), SURFACEAREA_MYELIN{i});              % nlamxnintn matrix
    LEAKMY_INODE(i, :) = 1./(sum(1./LEAKmy_INODE{i}, 1));                                                                         % 1xnitnn vec
end
LEAKI_INODE = unitsabs(par.intn.elec.pas.cond.units)*par.intn.elec.pas.cond.value.vec.*SURFACEAREA_INODE;                % vec

leak = zeros(tns, 2);
leak(nodes, 1) = reshape(LEAK_NODE', [], 1);
leak(intn, :) = [reshape(LEAKI_INODE', nis*nintn, 1), reshape(LEAKMY_INODE', nis*nintn, 1)];
clear LEAK_NODE LEAKmy_INODE LEAKMY_INODE LEAKI_INODE



% leak reversal potential
erev = unitsabs(par.elec.pas.erev.units)*par.elec.pas.erev.value.vec;






% longitudinal resistances
RaxNODE = unitsabs(par.node.elec.pas.axres.units)*par.node.elec.pas.axres.value.vec.*LENGTH_NODE./AREA_NODE;
RaxINODE = unitsabs(par.node.elec.pas.axres.units)*par.intn.elec.pas.axres.value.vec.*LENGTH_INODE./AREA_INODE;   % vec
RpINODE = unitsabs(par.myel.elec.pas.axres.units)*par.myel.elec.pas.axres.value.vec.*LENGTH_INODE./AREA_PERI;   % vec

Raxial = zeros(tnf,2);
Raxial(node2intn,2)=RpINODE(:,1)/2;
Raxial(intn2node,2)=RpINODE(:,nis)/2;
Raxial(intn2intn,2)=reshape(((RpINODE(:,2:end)+RpINODE(:,1:end-1))/2)',nintn*(nis-1),1);
Raxial(node2node,2)=0;
Raxial(node2intn,1)=RaxNODE(1:end-1, end)/2+RaxINODE(:,1)/2;
Raxial(intn2node,1)=RaxNODE(2:end, 1)/2+RaxINODE(:,end)/2;
Raxial(intn2intn,1)=reshape(((RaxINODE(:,2:end)+RaxINODE(:,1:end-1))/2)',nintn*(nis-1),1);
Raxial(node2node,1)=reshape(((RaxNODE(:,2:end)+RaxNODE(:,1:end-1))/2)',nnodes*(nns-1),1);










%% active conductances
Vbounds = [-200, 400];
Vstep = 0.01;
V = (Vbounds(1):Vstep:Vbounds(2))';

rateTableSize = length(V);
if ~isempty(par.channels)
    rateTable = generateRateFunctionTables(par.channels, V, dt, par.sim.temp);
end


for i = length(par.channels): -1 : 1
    activeconductance{i} = unitsabs(par.channels(i).cond.units)*par.channels(i).cond.value.*surfaceAreaAxolemma(par.channels(i).location);
end
for i = length(par.channels): -1 : 1
    gateVariable(i) = gatesSteadyState(par.channels(i), vrest(par.channels(i).location), par.sim.temp);
end

activesum = zeros(size(par.active.segments));
activesum2 = zeros(size(par.active.segments));








% -------------------------- CREATE "A" MATRIX ------------------------- %
A = spalloc(2*tns,2*tns,4*tns+4*(tns-1));
V1 = [vrest; zeros(tns, 1)];
V2 = [zeros(tns+2, 1), [[0,0]; reshape(V1,tns,2);[0,0]], zeros(tns+2,1)];

Gaxialpad = [[0,0]; 1./(2*Raxial); [0,0]];
leakpad=[zeros(tns,1),leak];
cappad=[zeros(tns,1),cap];

Radial=sum(cat(3,cappad,cappad(:,[2 3 1])),3)/dt + sum(cat(3,leakpad,leakpad(:,[2 3 1])),3)/2;
Radial=Radial(:,1:end-1);

Longitudinal=Gaxialpad(1:end-1,:)+Gaxialpad(2:end,:);

Total=Radial+Longitudinal;
Total(nodes,2)=ones(length(nodes),1);

Total(Total==Inf)=1;

Longitudinal2=-1./(2*Raxial);
% added to program - required
Longitudinal2([node2intn, node2node, intn2node],2)=0;

Radial2 = -(cap(:,1)/dt+leak(:,1)/2);
Radial2(nodes) = 0;


d = [-tns, -1, 0, 1, tns];

B = [[Radial2;zeros(tns,1)], ...
    reshape([Longitudinal2;[0,0]],2*tns,1), ...
    reshape(Total,2*tns,1), ...
    reshape([[0,0];Longitudinal2],2*tns,1), ...
    [zeros(tns,1);Radial2]];

A = spdiags(B, d, 2*tns, 2*tns);



Radialpre=(cappad/dt-leakpad/2);
Radialpre(isnan(Radialpre))=0;
Leakpre=[leak(:,1).*erev,-leak(:,1).*erev];

notcoveredAll  = nodes;
notcoveredAndActive = sort(intersect(par.active.segments, nodes), 'ascend');
covered = sort(intersect(par.active.segments, intn), 'ascend');
notcoveredActiveIdx = [];
for i = 1 : length(notcoveredAndActive)
    notcoveredActiveIdx(i, 1) = find(par.active.segments == notcoveredAndActive(i));
end
coveredActiveIdx = [];
for i = 1 : length(covered)
    coveredActiveIdx(i, 1) = find(par.active.segments == covered(i));
end

if isempty(notcoveredActiveIdx)
    notcoveredAndActive = [];
    activeUpdateNotCoveredIdx = [];
    offsetNotCovered = [];
else
    activeUpdateNotCoveredIdx = sub2ind(size(A), notcoveredAndActive, notcoveredAndActive);
    offsetNotCovered = A(activeUpdateNotCoveredIdx);
end

if isempty(coveredActiveIdx)
    covered = [];
    activeUpdateCoveredInsideIdx = [];
    offsetCoveredInside = [];
    activeUpdateCoveredInside2Idx = [];
    offsetCoveredInside2 = [];
    activeUpdateCoveredOutsideIdx = [];
    activeUpdateCoveredOutside2Idx = [];
    offsetCoveredOutside = [];
    offsetCoveredOutside2 = [];
else
    activeUpdateCoveredInsideIdx = sub2ind(size(A),covered,covered);
    offsetCoveredInside = A(activeUpdateCoveredInsideIdx);
    activeUpdateCoveredInside2Idx = sub2ind(size(A),covered,covered+tns);
    offsetCoveredInside2 = A(activeUpdateCoveredInside2Idx);
    activeUpdateCoveredOutsideIdx = sub2ind(size(A), covered+tns, covered+tns);
    activeUpdateCoveredOutside2Idx = sub2ind(size(A), covered+tns, covered);
    offsetCoveredOutside = A(activeUpdateCoveredOutsideIdx);
    offsetCoveredOutside2 = A(activeUpdateCoveredOutside2Idx);
end




if options.measureEnergy
    INaSave = zeros(T,1);
    energyNode = str2double(options.energyNode);
    naChannelNos = find(par.node.elec.act.issodium);
    numberOfSodiumChannels = length(naChannelNos);
end


if options.saveWithinModel
	fwrite(fid, V2(1+options.segmentsToSave, 2:3), 'double');
end
if ~options.saveWithinModel
    Vsave = vrest .* repmat([1, 0], tns, 1, T);
end

timetaken = toc;
fprintf('finished loading parameters, took %.2f s\n', timetaken);
fprintf('running...'); tic;


% --------------------------------------------------- MAIN ------------------------------------------------------ %

for i = 1:T
    i
    
    rateTableIdx = round(1 + rateTableSize * ((V2(2:end-1, 2) - Vbounds(1))/(Vbounds(2) - Vbounds(1))));
    if V2(2, 2) > Vbounds(2)
            disp('Too much stimulation: reduce stimulation')
        break
    end
    
    for ii = 1 : length(par.channels)
        gateVariable(ii).inf = rateTable(ii).gamma(rateTableIdx(par.channels(ii).location), :) + rateTable(ii).delta(rateTableIdx(par.channels(ii).location), :).*gateVariable(ii).inf;
    end
    
    
    activesum(:) = 0;
    activesum2(:) = 0;
    for ii = 1 : length(par.channels)
        temp = activeconductance{ii} .* prod(bsxfun(@power, gateVariable(ii).inf, par.channels(ii).gates.numbereach), 2);
        activesum(par.channels(ii).locationInActive) = activesum(par.channels(ii).locationInActive) + temp/2;
        activesum2(par.channels(ii).locationInActive) = activesum2(par.channels(ii).locationInActive) + temp.*par.channels(ii).erev.value;
    end
    
    INa = 2*activesum.*V2(1+par.channels(1).location, 2) - activesum2;
    IL = leak(par.channels(1).location, 1).*(V2(1+par.channels(1).location, 2) - erev(par.channels(1).location));
    
    A(activeUpdateNotCoveredIdx) = offsetNotCovered + activesum(notcoveredActiveIdx);
    A(activeUpdateCoveredInsideIdx) = offsetCoveredInside + activesum(coveredActiveIdx);
    A(activeUpdateCoveredInside2Idx) = offsetCoveredInside2 - activesum(coveredActiveIdx);
    A(activeUpdateCoveredOutsideIdx) = offsetCoveredOutside + activesum(coveredActiveIdx);
    A(activeUpdateCoveredOutside2Idx) = offsetCoveredOutside2 - activesum(coveredActiveIdx);
    
    
    V1 = Istim(:, :, (i>DUR+1)+1)...
            +(V2(1:end-2,2:3)-V2(2:end-1,2:3)).*Gaxialpad(1:end-1,:)...
            -(V2(2:end-1,2:3)-V2(3:end,2:3)).*Gaxialpad(2:end,:)...
            +Radialpre(:,1:2).*(diff(V2(2:end-1,1:3),1,2))...
            -Radialpre(:,2:3).*(diff(V2(2:end-1,2:4),1,2))...
            +Leakpre;
    
    
    V1(notcoveredAndActive) = V1(notcoveredAndActive)-activesum(notcoveredActiveIdx).*V2((tns+2) + (notcoveredAndActive+1))+activesum2(notcoveredActiveIdx);
    V1(notcoveredAll+tns) = 0;
     
    V1(covered) = V1(covered) - activesum(coveredActiveIdx).*(V2((tns+2) + (covered+1)) - V2(2*(tns+2) + (covered+1))) + activesum2(coveredActiveIdx);
    V1(covered+tns) = V1(covered+tns) + activesum(coveredActiveIdx).*(V2((tns+2) + (covered+1)) - V2(2*(tns+2) + (covered+1))) - activesum2(coveredActiveIdx);
%     
    V1 = V1(:); % reshape(V1,2*tns,1);
    
    % note: padarray is slower than the following
    V2 = [zeros(tns+2,1), [[0,0];reshape(A\V1,tns,2);[0,0]], zeros(tns+2,1)];
    
    if options.saveWithinModel
        fwrite(fid, V2(1+options.segmentsToSave, 2:3), 'double');
%         fwrite(fidINa, INa, 'double');
%         fwrite(fidIL, IL, 'double');
    else
        Vsave(:, :, i+1) = V2(2:end-1, 2:3);
    end
    
end



if options.saveWithinModel
    fclose(fid);
%     fclose(fidINa);
%     fclose(fidIL);
    Vsave = NaN;
end



% --------------------- END OF MAIN ---------------------- %

function out=constructshape(node,internodes,nis,tns)

temp = [node(1,1:end-1);repmat(internodes(1,:),nis,1)];
out(:,1)=[reshape(temp,tns-1,1);node(1,end)];
temp = [node(2,1:end-1);repmat(internodes(2,:),nis,1)];
out(:,2)=[reshape(temp,tns-1,1);node(2,end)];


