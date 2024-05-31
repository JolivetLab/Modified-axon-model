function g = gatesSteadyState(k, V, temperature)

if isempty(V)
    for channelIdx = 1 : length(k)
        g(channelIdx).inf = [];
    end
    return
end


rateFcn = generateRateFunctions(k);

for channelIdx = 1 : length(k)
    g(channelIdx).inf = nan([length(V), k(channelIdx).gates.number]);
    for j = 1 : k(channelIdx).gates.number
        a = rateFcn(channelIdx).alpha{j}(V, temperature);
        b = rateFcn(channelIdx).beta{j}(V, temperature);
        
        a = limitCases(V, a, rateFcn(channelIdx).alpha{j}, temperature);
        b = limitCases(V, b, rateFcn(channelIdx).beta{j}, temperature);
        
        g(channelIdx).inf(:, j) = a./(a + b);
    end
end