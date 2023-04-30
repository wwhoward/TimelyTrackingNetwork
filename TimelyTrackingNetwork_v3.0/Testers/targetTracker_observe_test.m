% Ensure targetTracker v3.0 performs as expected
% W.W.Howard, Wireless@VT, {wwhoward}@vt.edu

clc;clear;close all

myTargets = targetModel(10); 

myTracker = {targetTracker(1, 'Coverage', 2*pi)}; % coverage 2pi to guarentee total region coverage

t_step = 0.5; 
for t = 0:t_step:50
    t
    myTargets.update(t_step); 
    myTracker{1}.observe(myTargets, t); 
end

myTracker{1}.plotTrack(myTracker{1}.Targets{1}, myTargets.Targets{1})