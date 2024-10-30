function [spikes, t] = poissonSpikeGen(fr, t_step, n_neuron, sampling_rate)
%%% Inputs:
%%% fr: fire rate
%%% s_step: number of time pts
%%% n_neuron: number of neurons
%%% sampling_rate: temporal sampling rate, used to convert t_step to real time
%%% Outputs:
%%%     spikes: the spikes matrix
%%%     t: time
%%%%%%
%%% Codes taken from tutorial at https://praneethnamburi.com/author/praneethnamburi/

    dt = 1 / sampling_rate;
    tLength = t_step / sampling_rate;    
    spikes = rand(n_neuron, t_step) < fr*dt;
    t = 0:dt:tLength-dt;
end