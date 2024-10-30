function kernel = expDecay(tDecay, samplingRate)
%%% Get a exponential decay kernel
    tKernel = -8*tDecay : 1/samplingRate : 8*tDecay;
    kernelSize = size(tKernel, 2);
    kernel = exp(-tKernel/tDecay);
    kernel(1:floor(kernelSize/2)) = 0;
end