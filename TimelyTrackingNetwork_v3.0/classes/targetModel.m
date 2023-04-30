classdef targetModel < handle
    %targetModel simulates a group of multiple targets with shared properties
    %   Detailed explanation goes here
    
    properties % note only global properties, targets have their own as well
        % Cell array to hold target objects
        Targets

        % How many targets in each stage of life
        nActiveTargets      % Current active targets
        nRetiredTargets     % Current retired targets
        nTargets            % Current total targets

        % Desired mean # of targets
        saturation = 20

        % Spawn and retire probabiliities per second
        % Configured to hover around saturation
        pRetire = 0.01     % Influences target track duration
        pSpawn              % Initialized in constructor

        % Motion modeling
        motionModel         % Selects a particular set of motion models for all targets
        max_vel = 10        % Max velocity
        motionTypes = {'lin', 'turn', 'accel'}
        pLin = 0.9% This is hacky, but probability that motion stays linear. 

        % Regions
        spawnRegion = polyshape([100, 900, 900, 100], [100, 100, 900, 900])
        aliveRegion = polyshape([0, 1000, 1000, 0], [0, 0, 1000, 1000])
        
        % Keep track of the current time
        time
    end
    
    methods(Access=public)
        function obj = targetModel(initialTargets, varargin)
            %TARGETMODEL Handles target motion modeling
            
            obj.nTargets = initialTargets; % initial total targets
            obj.time = 0; % zero-indexing time

            obj.motionModel = 'lin_turns'; % {'lin', 'lin_turns', 'lin_turns_accel'}

            if any(strcmp(varargin, 'pLin'))
                idx = find(strcmp(varargin, 'pLin')); 
                obj.pLin = varargin{idx+1}; 
            end

            % Populate initial targets
            for i = 1:initialTargets
                obj.Targets{i} = obj.newTarget(obj.time); 
            end
            obj.update(0);

            obj.pSpawn = obj.pRetire*obj.saturation; 
        end % end constructor

        function [obj, flags] = update(obj, t)
            obj.time = obj.time + t; 
            flags = zeros(1, obj.nTargets); 

            % Update existing targets
            for i = 1:obj.nTargets
                if obj.Targets{i}.isActive
                    [obj.Targets{i}, flags(i)] = obj.updateTarget(obj.Targets{i}, t); 
                end % end if
            end % end for

            % Check for retirees to update internal tracking
            activeTargets = 0; 
            for i = 1:obj.nTargets
                if obj.Targets{i}.isActive
                    activeTargets = activeTargets + 1; 
                end % end if
            end % end for
            obj.nActiveTargets = activeTargets; 

            % Generate new targets according to Poisson
            newTargs = poissrnd(t*obj.pSpawn); 
            for i = 1:newTargs
                obj.Targets{end+1} = obj.newTarget(obj.time); 
                obj.nActiveTargets = obj.nActiveTargets + 1; 
                flags(end+1) = 1; 
            end
            obj.nTargets = size(obj.Targets, 2); 
            
            % Push flag updates back to each target
            for i = 1:obj.nTargets
                obj.Targets{i}.FlagHistory(end+1) = flags(i); 
            end

            % Calc nRetired
            obj.nRetiredTargets = obj.nTargets - obj.nActiveTargets; 
        end % end update      

        % Getters
        function [state, idx, tracks] = getState(obj, activityType)
            % Returns target states by career status
            idx = []; 
            tracks = {}; 
            switch activityType
                case "all"
                    state = zeros(4, obj.nTargets); 
                    for i = 1:obj.nTargets
                        state(:,i) = obj.Targets{i}.state; 
                        tracks{end+1} = obj.Targets{i}.Track; 
                        idx(end+1) = i; 
                    end % end for

                case "active"
                    state = zeros(4, obj.nActiveTargets); 
                    j=1; 
                    for i = 1:obj.nTargets
                        if obj.Targets{i}.isActive
                            state(:,j) = obj.Targets{i}.state; 
                            tracks{end+1} = obj.Targets{i}.Track; 
                            j=j+1; 
                            idx(end+1) = i; 
                        end % end if
                    end % end for

                case "retired"
                    state = zeros(4, obj.nRetiredTargets); 
                    j=1; 
                    for i = 1:obj.nTargets
                        if ~obj.Targets{i}.isActive
                            state(:,j) = obj.Targets{i}.state; 
                            tracks{end+1} = obj.Targets{i}.Track; 
                            j=j+1; 
                            idx(end+1) = i; 
                        end % end if
                    end % end for
            end % end switch
        end % end getState

        function [mean_age, ages] = getAge(obj)
            ages = []; 
            for i = 1:obj.nTargets
                ages(end+1) = obj.Targets{i}.Age; 
            end % end for

            mean_age = mean(ages); 
        end % end getAge


        % Setters
        function [] = set_pRetire(obj, pRetire)
            obj.pRetire = pRetire; 
            obj.pSpawn = obj.pRetire * obj.saturation; 
        end
    end % end public methods

    methods(Access=protected)
        % For updating
        function [targ, flag] = updateTarget(obj, targ, t)
            flag = 0; % unless it isn't
            targ.Age = targ.Age + t; 

            % First things: Update motion type
            current_idx = find(strcmp(obj.motionTypes, targ.Motion)); 
            n_models = length(targ.Transition(current_idx, :)); 
            trans = targ.Transition(current_idx, :); 
            trans(setdiff(1:n_models, current_idx)) = t*trans(setdiff(1:n_models, current_idx)); 
            trans(current_idx) = 1-sum(trans(setdiff(1:n_models, current_idx))); 
            new_idx = randsample(1:n_models, 1, true, trans); 
            if new_idx ~= current_idx
                flag = 1; % change in motion model raises flag
                targ.Motion = obj.motionTypes{new_idx}; 
                
                if current_idx == 1 && new_idx == 2
                    % Transition 1: const vel to const turn
                    % Update turn rate
                    targ.turn_rate = -1 + 2*rand; %Somewhere from -1:1
                elseif current_idx == 2 && new_idx == 1
                    % Transition 2: const turn to const vel
                    % Update vel
                    targ.state(2:2:end) = 0.1*randn*targ.state(2:2:end) + targ.state(2:2:end); % Randomly scale current velocity ~90:110%

                % TODO: Fill in the rest of the transitions for accel model

                end % end if
            end % end if

            % Second things: Update state
            switch targ.Motion
                case 'lin'
                    % Simply prop current vel

                    % Determine new position 
                    targ.state(1:2:end) = targ.state(1:2:end) + t*targ.state(2:2:end); 
                case 'turn'
                    % Get heading from vel, update pos, update vel

                    % Determine current heading
                    heading = atan2(targ.state(4), targ.state(2)); 

                    % Update heading according to turn rate
                    heading = heading + t*targ.turn_rate/1;

                    % Determine new velocity
                    speed = sqrt(targ.state(2)^2 + targ.state(4)^2); 
                    targ.state(2:2:end) = speed * [cos(heading), sin(heading)]; 

                    % Determine new position
                    targ.state(1:2:end) = targ.state(1:2:end) + t*targ.state(2:2:end); 
                case 'accel'
                    % TODO
            end % end switch

            % Third things: 
            % Can we retire yet? 
            if isinterior(obj.aliveRegion, targ.state([1,3]))
                if isinterior(obj.spawnRegion, targ.state([1,3]))
                    tmp_pRetire = obj.pRetire; 
                else
                    tmp_pRetire = 10*obj.pRetire; 
                end
            else
                tmp_pRetire = 100*obj.pRetire; 
            end
            if rand < t*tmp_pRetire
                targ.isActive = 0; 
                flag = 1; 
            end % end if 

            targ.Track(:, end+1) = targ.state; 
        end % end updateTarget


        
        % For initializing
        function targ = newTarget(obj, initTime)
            targ = struct(); 
            targ.initialTime = initTime; 
            targ.isActive = 1; 
            targ.Age = 0; 

            % Initialize motion model
            [targ.Motion, targ.Transition] = obj.getMotionModel; 

            % Get initial position and velocity (these change according to
            % motion model)
            pos = 800*rand(1,2)+100; % give some buffer so they can move out of spawn
            theta = 2*pi*rand; 
            speed = obj.max_vel * rand; 
            vel = speed * ([cos(theta), -sin(theta); sin(theta), cos(theta)] * ones(2,1)).'; 
            targ.turn_rate = -1 + 2*rand;

            targ.state = zeros(1,4); % this updates in each time step to always be current
            targ.state([1,3]) = pos; 
            targ.state([2,4]) = vel; 

            targ.FlagHistory = []; % times at which MM changes
            
            targ.Track(:,1) = targ.state; % State history
        end % end newTarget

        function [MM, Transition] = getMotionModel(obj)
            % Decides what type of motion is currently happening (not
            % specific values)
            switch obj.motionModel
                case 'lin' % const vel with random changes in dir & vel
                    % TODO
                case 'lin_turns' % const vel and const turn
                    % p1 = 0.8 + 0.2 * rand; % flat from 1/2 to 1
                    p1 = obj.pLin; 
                    p2 = 0.25 + 0.5 * rand; % flat from 1/4 to 3/4
                    % Transition = [p1, 1-p1; 1-p2, p2]; 
                    % Transition = [p1, 1-p1; 1-p1, p1]; 
                    Transition = [p1, 1-p1; 0.75, 0.25]; 
                    if rand > 0.5 % initialize to const vel
                        MM = 'lin';                         
                    else
                        MM = 'turn';
                    end % end if
                case 'lin_turns_accel' % const vel, const turn, const accel
                    % TODO
            end % end switch
        end % end getMotionModel
    end % end protected methods

    methods(Static)
        % Plotty boys
        function plotTrack(targ)
            track = targ.Track; 
            xy = track([1,3],:); 
            
            figure; 
            plot(xy(1,:), xy(2,:))
            legend('Track (Duration = ' + string(targ.Age) + ')', 'interpreter', 'latex', 'fontsize', 12)
            title('Ground Track', 'interpreter', 'latex', 'fontsize', 16)
            xlabel('x, meters', 'interpreter', 'latex', 'fontsize', 12)
            ylabel('y, meters', 'interpreter', 'latex', 'fontsize', 12)
        end % end plotTrack

        function plotnTargets(nTargets)
            figure; 
            plot(nTargets); 
            hold on
            yline(mean(nTargets)); 
            xlabel('Number of Active Targets', 'interpreter', 'latex', 'fontsize', 12)
        end % end plotnTargets
    end % end static methods
end % end class

