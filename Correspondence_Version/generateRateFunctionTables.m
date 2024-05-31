function rateTable = generateRateFunctionTables(k, V, dt, T)


rateFcn = generateRateFunctions(k);
for i = 1 : length(k)
    for j = 1 : k(i).gates.number
        a = rateFcn(i).alpha{j}(V, T);
        b = rateFcn(i).beta{j}(V, T);
        
        a = limitCases(V, a);
        b = limitCases(V, b);
        
        rateTable(i).gamma(:, j) = a./ (1/dt+(1/2)*(a + b));
        rateTable(i).delta(:, j) = (1/dt-(1/2)*(a + b))./(1/dt+(1/2)*(a + b));
    end
end


if isempty(k)
    rateTable(1).gamma = zeros(length(V), 1);
    rateTable(1).delta = zeros(length(V), 1);
end



function Y = limitCases(X, Y)
I = isinf(Y)|isnan(Y);
if any(I)
    Y(I) = interp1(X(~I), Y(~I), X(I));
end