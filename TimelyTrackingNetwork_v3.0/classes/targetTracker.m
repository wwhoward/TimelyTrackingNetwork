classdef targetTracker < handle
    % TARGETTRACKER Written for TTT Journal by W.W.Howard in Spring 2023
    %   Functions as a cognitive radar node. Works with 'fusionCenter' and
    %   'targetModel'. 
    % Version information: 
    % v3.0
    % 
    % Primarily for the study of Age of Information metrics as applied to
    % Cognitive Radar Networks. 
    % 
    % Also built with mode control in mind... 
    % 
    % Contact: {wwhoward}@vt.edu Wireless@VT
    
    properties
        nTargets = 0        % Estimate of how many targets are observable

        Targets = {}        % Set of all targets ever observed
        ActiveTargets = []  % Set of targets recently observed

        Age = []            % Age for each track. Stops updating when isActive=0. 
        ageThreshold = 5    % TODO
        innovationThreshold = 5 % TODO
        masterUpdateFlag = 0
        masterUpdateAge = 0
        masterUpdateAgeMax = 3

        ObservableRadius    % Determine how target is observed
        regionCoverage = 0.1 % Percent of region covered by this node
        NodePosition        % Self explanitory

        ID                  % Unique to this node

        delta_bar           % Constant for Bellman optimal policy

        time = 0

        updateChannels
        updateMethod = 'Centralized' 
        updateRate = 0.1; 
    end
    
    methods(Access=public)
        function obj = targetTracker(ID, varargin)
            %TARGETTRACKER Tracks targetModel targets, works with fusionCenter
            % Reqired inputs: 
                % ID: unique ID for this node
            % Optional inputs: 
                % 'Coverage': value in [0,1] denoting percentage of region covered by this node 
                % Note due to random positioning, set 'Coverage' to 2pi to guarentee total region coverage
            
            obj.ID = ID; 

            % Parse optional inputs
            % Valid optional inputs: 'Coverage', 'updateRate'
            
            if ~isempty(varargin)
                if any(strcmp(varargin, 'Coverage')) % covered region
                    idx = find(strcmp(varargin, 'Coverage')); 
                    obj.regionCoverage = varargin{idx+1}; 
                end

                if any(strcmp(varargin, 'updateRate')) % node update rate
                    idx = find(strcmp(varargin, 'updateRate')); 
                    obj.updateRate = varargin{idx+1}; 
                end
            end % end if
            
            % Define observable region
            obj.NodePosition = 1000*rand(2, 1); % Somewhere in one square km
            
            obj.ObservableRadius = sqrt(obj.regionCoverage * 1000^2/pi); 
        end

        function [flags, Targets] = observe(obj, targetModel, time)
            % targetModel: targets needing observing
            % time: current time

            obj.time = time; 
            flags = zeros(1, obj.nTargets); 

            [state, idx] = targetModel.getState('active'); 

            % Iterate through all possible new target IDs 
            % Target ID is constant so any new targets will have ID > length(obj.Targets)
            % Also, all target states & idx are provided, so we need to determine if each is observable
            % TODO: Replace this clunky struct with a class
            for i = length(obj.Targets)+1:max(idx)
                % if this i is still active, instance target for it
                % Only way i is not active is if it retired before we observed it
                obj.Targets{i} = {}; 
                obj.Targets{i}.Filter = obj.initializeFilter(state(:,idx==i)); 
                obj.Targets{i}.Track = []; % Raw observations
                obj.Targets{i}.FilteredTrack = []; % Filtered observations
                obj.Targets{i}.isActive = 1; 
                obj.Targets{i}.FlagHistory = []; % "something interesting happened"
                obj.Targets{i}.Age = 0; 
                obj.Targets{i}.Penalty = obj.time; 
                obj.Targets{i}.V = 0; % Last updated at t=0
                obj.Targets{i}.State = 0; 
                obj.Targets{i}.StateHat = -1; 
                obj.Targets{i}.ModelProb = [0.5, 0.5]; 
                obj.Targets{i}.Transitions = zeros(2,2); 
                obj.Targets{i}.TransitionProb = 0.5*ones(2,2); 
                obj.Targets{i}.IsNew = 1; % lets FC know we're new
                obj.Targets{i}.UpdateTimes = time; 

                flags(i) = 1; % If there's a new target, flag it! 
                obj.ActiveTargets(end+1) = i;                 
            end % end for

            % Now, iterate through the targets and update each one if it's observable
            transition_inactive = []; 
            transition_active = []; 
            for i = 1:length(idx)
                % Check for observability
                dist = norm(obj.NodePosition - state([1,3],i)); 
                if dist > obj.ObservableRadius % not in range
                    % Check for active -> inactive transition
                    if ~isempty(intersect(obj.ActiveTargets, idx(i))) % transition case
                        transition_inactive(end+1) = idx(i); 
                        obj.Targets{idx(i)}.isActive = 0; % Set as inactive so it won't be in updates
                        obj.ActiveTargets = setdiff(obj.ActiveTargets, idx(i)); % Remove from active targets
                    end % end if                    
                else    
                    % Check for inactive -> active transition
                    if isempty(intersect(obj.ActiveTargets, idx(i))) % transition case
                        transition_active(end+1) = idx(i); 
                        obj.ActiveTargets = union(obj.ActiveTargets, idx(i)); % Add to active targets if in range
                        obj.Targets{idx(i)}.isActive = 1;
                    end % end if
    
                    % Some notes on target tracking: 
                    % obj.Targets{}.Filter is an Interacting Motion Model filter
                    % It evaluates between three different possible motion models 
                    %   cv, ct, ca
                    % 
                    % obj.Targets{}.Track stores the raw observed [x,y]
                    % pstates is the [x vx y vy z vz] prediction
                    % cstates is the [x vx y vy z vz] correction
                    % obj.Targets{}.FilteredTrack stores the corrected [x,y]
                    % The filter is 3D because Old Will couldn't 2D it
    
                    % The variance for a track is based on the true range
                    % We know that it's only detected if 0<dist<100 so...
                    variance = (dist/obj.ObservableRadius)*10; 
    
                    % Since targets may enter and exit the region ... 
                    td = time - obj.Targets{idx(i)}.UpdateTimes(end); 
                    obj.Targets{idx(i)}.UpdateTimes(end+1) = time; 
                    
    
                    obj.Targets{idx(i)}.Track(:,end+1) = variance*randn(2,1) + state([1,3],i); 
    
                    pstates = predict(obj.Targets{idx(i)}.Filter, td); 
                    cstates = correct(obj.Targets{idx(i)}.Filter, [obj.Targets{idx(i)}.Track(:,end); 0]); 
                    obj.Targets{idx(i)}.FilteredTrack(:,end+1) = cstates([1,2,3,4]); 
    
                    % Update transition probabilities
                    tmpModelProbs = obj.Targets{idx(i)}.Filter.ModelProbabilities([1,3]); 
                    obj.Targets{idx(i)}.ModelProb(end+1, :) = tmpModelProbs ./ sum(tmpModelProbs); % renormalize
                    [~, states] = max(obj.Targets{idx(i)}.ModelProb(end-1:end,:), [], 2); 
                    tmp = zeros(2,2); 
                    tmp(states(1), states(2)) = 1; 
                    obj.Targets{idx(i)}.Transitions = obj.Targets{idx(i)}.Transitions + tmp; 
                    if ~any(sum(obj.Targets{idx(i)}.Transitions, 2)==0)
                        obj.Targets{idx(i)}.TransitionProb = obj.Targets{idx(i)}.Transitions./sum(obj.Targets{idx(i)}.Transitions, 2); 
                    end % end if
                    obj.Targets{idx(i)}.State = states(2); 
                    
                    % Flag high innovation (old method, here for legacy
                    % reasons. 
                    % TODO update this probably 
                    if vecnorm(pstates([1,3]) - obj.Targets{idx(i)}.Track(:,end),2,1)>obj.innovationThreshold; 
                        flags(idx(i)) = 1; 
                    end % end if           
                end % end ifelse
                flags(transition_active) = 1; % Flag active transitions
                flags(transition_inactive) = 1; % Flag inactive transitions
            end % end for
            Targets = obj.Targets; % Optionally output all target tracks
        end % end observe
        
        function [] = update(obj)
            obj.updateChannels.newUpdate(obj.time, obj.Targets, obj.ActiveTargets)
        end % end update

        % Setters
        function [] = setUpdateMethod(obj, method)
            % Determines when targetTracker posts updates to updateChannel
            switch method
                case 'centralized'
                    obj.updateMethod = @obj.sendCentralizedUpdates; 
                case 'distributed_random'
                    obj.updateMethod = @obj.sendRandomUpdates; 
                case 'distributed'
                    obj.updateMethod = @obj.sendAoiiUpdates; 
            end % end switch
        end % end setUpdateMethod
        

        % Getters
        function [state, idx] = getState(obj)
            % getState: returns current Kalman state for all targets
            % recently active
            % "recently active" -> ismember(obj.ActiveTargets)

            state = zeros(4, length(obj.ActiveTargets)); 
            idx = obj.ActiveTargets; 

            for i = 1:length(obj.ActiveTargets)
                pstates = predict(obj.Targets{obj.ActiveTargets(i)}.Filter, 0); 
                state(:,i) = pstates; 
            end % end for
        end % end getState

        function [targets, idx] = getTargets(obj, status)
            switch status 
                case 'all'
                    targets = obj.Targets{obj.ActiveTargets}; 
                    idx = obj.ActiveTargets; 
                case 'active'
                    targets = obj.Targets; 
                    idx = 1:length(obj.Targets); 
            end
        end

    end % end public methods

    methods(Access=private)
        function [] = sendCentralizedUpdates(obj)
            % Do nothing; FC prompts updates when necessary
        end % end sendCentralizedUpdates

        function [] = sendRandomUpdates(obj)
            % Randomly sends an update with a mean of updateRate
            if rand < obj.updateRate
                obj.update(); 
            end % end if
        end

        function [] = sendAoiiUpdates(obj)

        end


    end % end private methods

    methods(Static)
        % Internal magic (don't let the smoke out) 
        function filter = initializeFilter(state)
            state = [state([1,3]); 0]; 
            filter = trackingIMM(); 
            [~] = predict(filter, 0); 
            [~] = correct(filter, state); 
        end % end initializeFilter
        
        % Plotters
        function plotTrack(estTarget, Target)
            est_track = estTarget.Track; 
            filt_track = estTarget.FilteredTrack; 
            track = Target.Track; 
            
            figure; 
            plot(track(1,:), track(3,:)); 
            hold on
            plot(est_track(1,:), est_track(2,:), '.'); 
            plot(filt_track(1,:), filt_track(3,:), 'linewidth', 2); 
            % plot(track(1,Target.FlagHistory(1:size(track,2))==1), track(3,Target.FlagHistory(1:size(track,2))==1), 'x'); 
            % plot(filt_track(1,estTarget.FlagHistory(1:size(est_track,2))==1), filt_track(3,estTarget.FlagHistory(1:size(est_track,2))==1), 'o')
            legend('Track, age=' + string(Target.Age), 'Estimated Track', 'Filtered Track', 'Events', 'Flags', 'interpreter', 'latex', 'fontsize', 12); 
        end % end plotTrack
    end % end static methods
end

