
% [ sttc , DT ] = maksttc_cutts ( win , dt , A , B )
% 
% A MEX adaptation of the spike_time_tiling_coefficient.c function provided
% by Cutts and Eglen at https://github. com/CCutts/
% Detecting_pairwise_correlations_in_spike_trains. This uses an O( n ^ 2 )
% algorithm. The MEX wrapper will evaluate sttc between spike trains A and 
% B within analysis window win at each delta-t value from 0 to dt, rounded
% up to the nearest millisecond. This adds another level of nested looping,
% so the overal algorithm is O( n ^ 3 ). A and B must be double floating
% point vectors of spike times in chronological order, and dt must be a
% scalar double with delta-t in seconds. win is a two-element double of the
% start and end of the window in seconds. Returns a double vector. Optional
% output DT is a list of delta-t times in register with sttc, also a double
% vector.
% 
% It is added to MAK mainly for validation of maksttc.
% 
% Reference: 
% 
%   Cutts CS, Eglen SJ. 2014. Detecting Pairwise Correlations in Spike
%     Trains: An Objective Comparison of Methods and Application to the
%     Study of Retinal Waves. J Neurosc, 34(43):14288-14303.
% 
% Adapted by Jackson Smith - March 2018 - DPAG , University of Oxford
