function erev = updateleakerev(par)

vrest = unitsabs(par.elec.pas.vrest.units) * par.elec.pas.vrest.value.vec;
Iactive = zeros(par.geo.totalNumberSegments, 1);
for i = 1:length(par.channels)
    
    if isempty(par.channels(i).location), continue, end
    
    gateVariable(i) = gatesSteadyState(...
        par.channels(i), ...
        vrest(par.channels(i).location), ...
        par.sim.temp);
    
    Iactive(par.channels(i).location) = Iactive(par.channels(i).location) + ...
        currentInstantaneous(par.channels(i), vrest(par.channels(i).location), gateVariable(i));
    
end

erev = vrest + Iactive ./ (unitsabs(par.elec.pas.cond.units)*par.elec.pas.cond.value.vec);
