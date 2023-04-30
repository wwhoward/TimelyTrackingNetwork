classdef updateChannel < handle
    %UPDATECHANNEL Written for TTT Journal by W.W.Howard in Spring 2023
    %   Allows information exchange between targetTracker and fusionCenter
    
    properties
        updateTimes = zeros(1)
        updates = {}
    end
    
    methods
        function obj = updateChannel()
            %UPDATECHANNEL Provides updates from targetTracker to fusionCenter
        end

        function [] = newUpdate(obj, time, tracks, active_targets)
            % newUpdate pushes time and tracks to updateChannel
            obj.updateTimes(end+1) = time; 
            obj.updates{end+1} = {tracks, active_targets}; 
        end
    end
end

