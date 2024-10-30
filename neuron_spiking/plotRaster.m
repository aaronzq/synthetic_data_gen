function [] = plotRaster(spikes, t)
%%% Visualize the spikes with a raster plot
%%% Codes taken from tutorial at https://praneethnamburi.com/author/praneethnamburi/
hold all;
for trialCount = 1:size(spikes,1)
    spikePos = t(spikes(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikes, 1)+1]);