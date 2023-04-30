function [trackers] = genTrackers(nTrackers, varargin)
    % GENTRACKERS Convenience function to enable instancing of targetTrackers in parfor loops
    % Written for TTT Journal by W.W.Howard in Spring 2023
    % Version information: 
    % v3.0
    % Contact: {wwhoward}@vt.edu Wireless@VT
    
    trackers = {}; 
    for n = 1:nTrackers
        trackers{n} = targetTracker(n, varargin{:}); 
    end
end