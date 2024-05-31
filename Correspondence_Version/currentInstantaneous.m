function I = currentInstantaneous(k, V, M)

% V - voltage in mV
% I - current in uA/mm^2

I = zeros(size(V));

for i = 1 : length(k)
    g = prod(bsxfun(@power, M(i).inf, k(i).gates.numbereach), 2);
    I = I + unitsabs(k(i).cond.units)*k(i).cond.value .* g .* (V - unitsabs(k(i).erev.units)*k(i).erev.value);
end