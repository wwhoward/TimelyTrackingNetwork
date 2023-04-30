function [averageStats] = AverageStats(stats)
%AVERAGESTATS Written for TTT Journal by W.W.Howard in Spring 2023
% Version information: 
% For TimelyTrackingNetwork v3.0
% Contact: {wwhoward}@vt.edu

% Input: 
% Stats -> (nParam by nRep)

% Output: 
% averageStats -> stats, averaged over the nRep dimension

averageStats = dictionary(); 
nParams = size(stats, 1); 
nRep = size(stats, 2); 


% Get initial TimeSteps for comparison
baseTimeSteps = stats{1,1}{"TimeSteps"}; 
nSteps = length(baseTimeSteps); 

age = zeros(nParams, nRep, nSteps); 
peak_age = zeros(nParams, nRep, nSteps); 
err = zeros(nParams, nRep, nSteps);
rmse  = zeros(nParams, nRep, nSteps);
covered  = zeros(nParams, nRep, nSteps);
tracked  = zeros(nParams, nRep, nSteps);
total = zeros(nParams, nRep, nSteps); 

for p = 1:nParams
    for n = 1:nRep
        if any(stats{p, n}{"TimeSteps"} ~= baseTimeSteps)
            error("Can't synchronize time steps... ")
        end
        age(p, r, :) = stats{p, n}{"Age"}; 
        peak_age(p, r, :) = stats{p, n}{"PeakAge"}; 
        err(p, r, :) = stats{p, n}{"Error"}; 
        rmse(p, r, :) = stats{p, n}{"RMSE"}; 
        covered(p, r, :) = stats{p, n}{"nCoveredTargets"}; 
        tracked(p, r, :) = stats{p, n}{"nTrackedTargets"}; 
        total(p, r, :) = stats{p, n}{"nTotalTargets"}; 
    end % end for nRep
end % end for nParams

averageStats("TimeSteps") = baseTimeSteps; 

averageStats("Age") = mean(age, 2); 
averageStats("PeakAge") = mean(peak_age, 2); 
averageStats("Error") = mean(err, 2); 
averageStats("RMSE") = mean(rmse, 2); 
averageStats("nCoveredTargets") = mean(covered, 2); 
averageStats("nTrackedTargets") = mean(tracked, 2); 
averageStats("nTotalTargets") = mean(total, 2); 


% Notes
% Need age, peak age, error, rmse, covered t, tracked t, total t




end

