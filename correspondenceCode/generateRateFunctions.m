function rateFcn = generateRateFunctions(channels)

rateNames = {'alpha', 'beta'};
for rateIdx = 1 : length(rateNames)
    for channelIdx = 1 : length(channels)
        for gateIdx = 1 : channels(channelIdx).gates.number
            rateFcn(channelIdx).(rateNames{rateIdx}){gateIdx} = ...
                eval(sprintf('@(V, T)((%.2f^((T-%.1f)/10))*(%s))', ...
                channels(channelIdx).gates.(rateNames{rateIdx}).q10(gateIdx), ...
                channels(channelIdx).gates.temp, ...
                channels(channelIdx).gates.(rateNames{rateIdx}).equ{gateIdx}));
        end
    end
end