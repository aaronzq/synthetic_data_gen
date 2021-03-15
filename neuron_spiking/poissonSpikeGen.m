function [spikeMat, tVec] = poissonSpikeGen(fr, tStep, nNeuron, samplingRate)    
    dt = 1 / samplingRate;
    tLength = tStep / samplingRate;    
    spikeMat = rand(nNeuron, tStep) < fr*dt;
    tVec = 0:dt:tLength-dt;
end