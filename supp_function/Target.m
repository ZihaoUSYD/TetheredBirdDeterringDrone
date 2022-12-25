classdef Target
    % ===============================================================
    % Class - Target
    % This class contains all the properties of the target object
    %
    % Methods: Protected
    %           [obj] = initInterests(obj, world)
    %           [obj] = integrateVelocity(obj, cumulation, dt)
    %           [obj] = integrateHeading(obj, cumulation, dt)
    %           [obj] = integratePosition(obj, cumulation, dt)
    %           [obj] = updateInterests(obj, agents, world, cumulation)
    %           [obj] = moveTarget(obj, world, cumulation, dt)
    %          Private
    %           [obj] = expandArrays(obj, world, iter)
    %
    % Created by: Zihao Wang
    % Last updated: 19-07-2018
    % ===============================================================
    properties
        % states
        x                           % x position
        y                           % y position
        z                           % z position
        heading                     % heading
        velocity                    % speed
		vertical_velocity			% change in altitude (decoupled from lateral dynamics)
        flight                      % flying or landed
		energy						% energy reserve
        
        % performance
        sensor_range                % how far can they sense agents
        max_yawrate                 % maximum rate of heading change (rad/s)
        max_velocity                % maximum velocity
        
        % interests map
        inter_map                   % the target's interested location
        
        % figure handles
        figure_handle_uncertainty   % figure handle for uncertainty
        figure_handle_probability   % probability figure handle
        figure_handle_interests     % interests figure handle
        figure_handle_waypoint      % waypoint plot handle
        figure_handle_debug         % figure handle for debug
        interests_map_handle        % interests map handle
        
        % system dynamics
        velocity_setpoint           % velocity setpoint
		vertical_velocity_setpoint	% vertical velocity setpoint
        heading_setpoint            % heading setpoint
        x_setpoint                  % position x setpoint
        y_setpoint                  % position y setpoint
		z_setpoint					% position z setpoint
		z_flying					% reference z setpoint for birds flying
		descent_flag				% Boolean for tracking whether bird is undergoing descent or not
		
        % gains
        pos_p               % position error to velocity setpoint gain
		vertPos_p			% vertical position error to vertical velocity setpoint gain
        vel_p               % velocity error to acceleration setpoint gain
        yaw_p               % heading error to yaw speed setpoint gain
		vertVel_p			% vertical velocity error to vertical acceleration setpoint gain
		
        % cost map (experimental)
        cost
        % grid target distance (experimental)
        grid_target_dist
		
		% state property indicating whether target has left property or not
		in_property
		in_property_before
		d_outside
		
		% bird name and other metrics
		bird_name
		BMR
		m_bird
		P_met_burst
		P_met_long

    end
    
    methods
        %% constructor
        function obj = Target(flight,x,y,z,heading,velocity,sensor_range,max_velocity,max_yawrate, bird_name)
            if nargin > 0
                obj.flight      = flight;
                obj.x           = x;
                obj.y           = y;
                obj.z           = z;
                obj.heading     = heading;
                obj.velocity    = velocity;
				obj.vertical_velocity = 0;
                obj.sensor_range= sensor_range;
                obj.max_velocity= max_velocity;
                obj.max_yawrate = max_yawrate;
                obj.x_setpoint  = x;
                obj.y_setpoint  = y;
				obj.z_setpoint	= z;
				obj.z_flying	= z;
                obj.energy		= 1020; % J, scale between 0-10, assume that the bird is not at maximum energy
				obj.bird_name	= bird_name;
				
            % assign default values if no values are given
            else
                obj.flight      = true;
                obj.x           = 0;
                obj.y           = 0;
                obj.z           = 0;
				obj.vertical_velocity = 0;
                obj.heading     = 0;
                obj.velocity    = 10;
                obj.sensor_range= 400;
                obj.max_velocity= 10;
                obj.max_yawrate = deg2rad(20);
                obj.x_setpoint  = 0;
                obj.y_setpoint  = 0;
				obj.energy		= 0.9;
				obj.bird_name	= 'Starling';
            end
            
            % gains
            obj.pos_p = 0.12;
			obj.vertPos_p = 0.1;
            obj.vel_p = 0.65;
			obj.vertVel_p = 0.6;%0.3%0.6; % birds appear to have very quick response times/abilities to alter their vertical velocity?
            obj.yaw_p = 0.5;
			
			% assume target starts inside property boundary and that it's
			% the first time
			obj.in_property = true;
			obj.in_property_before = false;
			obj.descent_flag = false; % assume that the birds are not undergoing descent straight away 
			obj.d_outside = 10; % starting distance to leave property
			% Set up metrics/data specific to target bird type
			% Access bird catalogue
			[numData,~,cellData] = xlsread('birdCatalogueData.xlsx');
			
			for i = 2:length(cellData(:,1)) % skip header names
				
				% Only extract data for desired bird
				if cellData{i,1} == obj.bird_name
					% NOTE: numData exludes all non-numeric types, hence it is
					% smaller, account for this
					obj.BMR = numData(i-1,1);
					obj.m_bird = numData(i-1,2);
					obj.P_met_burst = numData(i-1,3);
					obj.P_met_long = numData(i-1,4);
				end
			end
			
        end
        
        %% protected methods
        % ========================
        % initialise interests map
        % ========================
        function [obj] = initInterests(obj, world)
            obj.inter_map           = ones(world.number_of_tiles_y, world.number_of_tiles_x).*world.inter0;
		end
        		  
        % ========================
        % update interests
        % ========================
        function [obj] = updateInterests(obj, agents, world, cumulation)
            % load the last global interests into a temporary variable
            data                = world.inter_global(:,:,cumulation-1);

            % drive interests value towards nominal
            data                = data.*world.inter_tau + world.inter_nom*(1-world.inter_tau);
            
            % extract coordinates to make code look cleaner
            target_x            = obj.x(cumulation-1);
            target_y            = obj.y(cumulation-1);

            for i = 1:size(agents,1)
                % birds do not react to fixed wing type UAV or ground
                % cameras or undeployed agents
                if agents(i).type == 1 || agents(i).type == 3 || agents(i).state == "HOME"
                    continue;
                end

                % extract coordinates to make code look cleaner
                agent_x             = agents(i).x(cumulation-1);
                agent_y             = agents(i).y(cumulation-1);

                % calculate distance between target and agent
                target_agent_dist   = euclideanDistance(agent_x,agent_y,target_x,target_y);

                % if agent is not in target's sensor range, continue
                if target_agent_dist > obj.sensor_range 
                    continue;
                end

                % calculate distance between all grids and agent
                grid_agent_dist     = euclideanDistance(agent_x,agent_y,world.grid_centre_x,world.grid_centre_y);

                % all tiles in the alert radius have a decreased interests proportional to distance
                intersection_alert  = grid_agent_dist < world.inter_alert_radius;
                % minimum of the existing value or new value
                data(intersection_alert) = min(data(intersection_alert), grid_agent_dist(intersection_alert)./world.inter_alert_radius);

                % the rest of the tiles in the nogo radius have 0 interests
                intersection_nogo   = grid_agent_dist < world.inter_nogo_radius;
                data(intersection_nogo) = world.inter_nogo;
            end % end agents for loop
			
			% Invert interests if the target no longer has sufficient
			% energy
% 			if obj.energy < 0
% 				data = 1./data;
% 			end
            % update target
            obj.inter_map(:,:,cumulation) = data;
        end
		
		% ========================
        % move target
        % ========================
        function [obj] = moveTarget(obj, world, cumulation, dt)
            % extract coordinates to make code look cleaner
            target_x                = obj.x(cumulation-1);
            target_y                = obj.y(cumulation-1);
            
            % calculate distance between target and grid centres
            obj.grid_target_dist        = euclideanDistance(target_x,target_y,world.grid_centre_x,world.grid_centre_y);
            
            % find all grids unreachable by the target using receding horizon
            unreachable_grids       = ~(obj.grid_target_dist < (obj.max_velocity*(world.RHC_steps)*dt));
            
			% find all grids not within vision radius (i.e. birds can't
			% just choose the global interest maximum if they have no
			% knowledge of it)
			
			
            % calculate relative heading change between target and grid centres
            grid_target_yaw         = wrapTo2Pi(atan2((world.grid_centre_y-target_y),(world.grid_centre_x-target_x)));
            
            % calculate cost to reach each grid centres
            % weights
            alpha                   = 0.85; % interests weight, was 0.85
            beta                    = 0.3; % distance weight, was 0.3
            gamma                   = 0.35; % heading weight
            % interests cost
            interests_cost          = alpha*(1-world.inter_global(:,:,cumulation-1));
            % distance cost
            distance_cost           = beta*obj.grid_target_dist./max(max(obj.grid_target_dist));
            % heading cost
            heading_difference      = abs(grid_target_yaw-obj.heading(cumulation-1));
            heading_difference(heading_difference>=pi)=2*pi-heading_difference(heading_difference>=pi);
            heading_cost            = gamma.*heading_difference./pi;
            % unreachable grids
            unreachable_cost        = 5*unreachable_grids;
            % total cost
            obj.cost                = unreachable_cost + ...
                                        interests_cost + ...
                                        distance_cost + ...
                                        heading_cost;
                                  
            % find the optimum waypoint
            [temp,waypoint_x]       = min(obj.cost);
            [~,waypoint_y]          = min(temp);
            waypoint_x              = waypoint_x(waypoint_y);
            
			% find first/any 'suitable' grid 
% 			interest_threshold		= 0.1;
% 			suitable_grids			= obj.cost <= interest_threshold;	
% 			reduced_suitable_indices  = [];	
% 			
% 			% Extracting the indices of all grids within property that
% 			% satisfy the minimum cost threshold (i.e. maximum interest)
% 			for i = 1:length(suitable_grids(1,:))
% 				for j = 1:length(suitable_grids(:,1))
% 					if suitable_grids(i,j) == 1
% 						reduced_suitable_indices = [reduced_suitable_indices; i, j];
% 					end
% 				end
% 			end
% 			
% 			% Sample random choice from 'suitable' options
% 			if ~isempty(reduced_suitable_indices) 
% 				randIndx = randi(length(reduced_suitable_indices(:,1)),1);
% 				rand_waypoint_x = reduced_suitable_indices(randIndx,1);
% 				rand_waypoint_y = reduced_suitable_indices(randIndx,2);
% 			else % if no suitable choices, use old waypoint
% 				rand_waypoint_x = waypoint_x;
% 				rand_waypoint_y = waypoint_y;
% 			end
% 			
% 			% Weighted random sample between the 'optimal' grid vs. any
% 			% suitable grid
% 			grid_options = ["Optimal", "Random"];
% 			grid_option_weights = [0.8, 0.2]; % preferentially choose the optimal waypoint 80% of the time - arbitrarily chosen atm
% 			grid_choice = randsample(grid_options,1,true,grid_option_weights);
% 						
% 			% Randomise the acceptable threshold to begin descent to crop
% 			% (note that this currently updating on EVERY time step, not
% 			% just to initiate the descent sequence)
			lower_dist = 10;
			upper_dist = 30;
			randDistanceThreshold = lower_dist + (upper_dist).*rand();
% 			horizontalDistanceToRandWaypoint = euclideanDistance(target_x,target_y,world.grid_centre_x(rand_waypoint_x,rand_waypoint_y),world.grid_centre_y(rand_waypoint_x,rand_waypoint_y));
% 			
% % 			if cumulation == 10
% % 				gmmy = 0;
% % 			end
% 			
% 			% Allow change to random waypoint (if presented from weighted random sample) if it is suitably far away from current waypoint (since birds can be lazy)		
% 			if grid_choice == "Random" && (horizontalDistanceToRandWaypoint < randDistanceThreshold)
% 				% change waypoint to random sample (else it's already the
% 				% optimal choice)
% 				waypoint_x = rand_waypoint_x;
% 				waypoint_y = rand_waypoint_y;
% 			end
			
            % assign the optimum waypoint to the position set point
            obj.x_setpoint(cumulation) = world.grid_centre_x(waypoint_x,waypoint_y);
            obj.y_setpoint(cumulation) = world.grid_centre_y(waypoint_x,waypoint_y);
			
			% If bird no longer in property change setpoint to arbitrary
			% point 10m perpendicular to property boundary

			if ~obj.in_property
				rot_angle = wrapTo2Pi(atan2(obj.y(cumulation-1)-world.world_centre(2),obj.x(cumulation-1)-world.world_centre(1)));
				
				if rot_angle <= pi()/4 || rot_angle > 7/4*pi()
					% right side
					obj.x_setpoint(cumulation) = world.size_x+20;
					obj.y_setpoint(cumulation) = world.world_centre(2);
				elseif rot_angle > pi()/4 && rot_angle <= 3/4*pi()
					% top side
					obj.x_setpoint(cumulation) = world.world_centre(1);
					obj.y_setpoint(cumulation) = world.size_y + 20;
				elseif rot_angle > 3/4*pi() && rot_angle <= 5/4*pi()
					% left side
					obj.x_setpoint(cumulation) = -20;
					obj.y_setpoint(cumulation) = world.world_centre(2);
				else
					% bottom side
					obj.x_setpoint(cumulation) = world.world_centre(1);
					obj.y_setpoint(cumulation) = -20;
					
				end
				
				obj.z_setpoint(cumulation) = obj.z_flying;
				
% 				tmp_heading = wrapTo2Pi(obj.heading(cumulation-1));
% 				if tmp_heading == 0 || tmp_heading == 2*pi()
% 					obj.x_setpoint(cumulation) = obj.x(cumulation-1)+obj.d_outside;
% 					obj.y_setpoint(cumulation) = obj.y(cumulation-1);
% 				elseif tmp_heading == deg2rad(90)
% 					obj.x_setpoint(cumulation) = obj.x(cumulation-1);
% 					obj.y_setpoint(cumulation) = obj.y(cumulation-1)+obj.d_outside;
% 				elseif tmp_heading == deg2rad(180)
% 					obj.x_setpoint(cumulation) = obj.x(cumulation-1)-obj.d_outside;
% 					obj.y_setpoint(cumulation) = obj.y(cumulation-1);
% 				elseif tmp_heading == deg2rad(270)
% 					obj.x_setpoint(cumulation) = obj.x(cumulation-1);
% 					obj.y_setpoint(cumulation) = obj.y(cumulation-1)-obj.d_outside;
% 				else
% 					% regular cases
% 					obj.y_setpoint(cumulation) = obj.y(cumulation-1) + obj.d_outside*sin(tmp_heading);
% 					obj.x_setpoint(cumulation) = obj.x(cumulation-1) + (obj.y(cumulation)-obj.y(cumulation-1))/min(tan(tmp_heading),1);
% 				end
% 				obj.z_setpoint(cumulation) = 15;
% 				% Reduce distance incrementally
% 				if obj.d_outside - obj.velocity(cumulation-1)*dt >= 0
% 					obj.d_outside = obj.d_outside - obj.velocity(cumulation-1)*dt;
% 				else
% 					obj.d_outside = 0;
% 				end
			end
			% If bird is within certain range of desired interest and horizontal velocity is therefore decreasing, begin
			% descent, otherwise try to fly at 'maximum default' altitude
			horizontalDistanceToWaypoint = euclideanDistance(target_x, target_y, obj.x_setpoint(cumulation), obj.y_setpoint(cumulation));
			
			if cumulation < 3 % neeed more time steps for checking velocity change
				velocityChange = 0;
			else
				velocityChange = (obj.velocity(cumulation-1) - obj.velocity(cumulation-2))/dt; % check whether bird flock is decelerating or accelerating (i.e. being chased)
			end

			% Descent to ground if below random threshold for first time OR
			% bird flock has already begun descending - only valid for both
			% cases while decelerating
			if (horizontalDistanceToWaypoint <= randDistanceThreshold || obj.descent_flag) && velocityChange <= 0
				% Target seeks crop on ground
				obj.z_setpoint(cumulation) = 1.5; % Hard-coded, simulates shrub/crop height above ground
				obj.descent_flag = true; % set flag to true, since we have begun descending
			else
				% Target is not within range of desired or interest or it
				% has now been scared from the current position and must
				% fly away --> seek flying height 
				obj.z_setpoint(cumulation) = obj.z_flying;
				obj.descent_flag = false; % reset flag to false
			end
				
	
            % integrate; CAN ONLY do this action with sufficient time
            % history: once cumulation >= 5 
			% Similarly, if current target has left property stop moving it
			if cumulation >= 2 %&& obj.in_property && ~obj.in_property_before
				obj                     = integratePosition(obj, cumulation, dt);
			end
		end
		
		
		% ========================
        % integrate position
        % ========================
        function [obj] = integratePosition(obj, cumulation, dt)
            		
			% distance error = setpoint - current position
            distance_error          = sqrt((obj.x_setpoint(cumulation)-obj.x(cumulation-1))^2+(obj.y_setpoint(cumulation)-obj.y(cumulation-1))^2);
			vert_distance_error		= obj.z_setpoint(cumulation)-obj.z(cumulation-1);
			
            % velocity setpoint = position gain * distance error
            % constrain velocity setpoint to the maximum
            obj.velocity_setpoint(cumulation) = min(obj.pos_p*distance_error,obj.max_velocity);
			obj.vertical_velocity_setpoint(cumulation) = min(obj.pos_p*vert_distance_error, obj.max_velocity);			
			
            % integrate velocity
            [obj]                   = integrateVelocity(obj, cumulation, dt, "Horizontal"); % Horizontal
			[obj]                   = integrateVelocity(obj, cumulation, dt, "Vertical"); % Vertical
			
            % heading setpoint = angle between position setpoint and current position in global coordinates
            obj.heading_setpoint    = wrapTo2Pi(atan2((obj.y_setpoint(cumulation)-obj.y(cumulation-1)),obj.x_setpoint(cumulation)-obj.x(cumulation-1)));
            % integrate heading
            [obj]                   = integrateHeading(obj, cumulation, dt);
            % integrate position
            obj.x(cumulation)       = obj.x(cumulation-1) + obj.velocity(cumulation-1)*cos(obj.heading(cumulation-1))*dt;
            obj.y(cumulation)       = obj.y(cumulation-1) + obj.velocity(cumulation-1)*sin(obj.heading(cumulation-1))*dt;
			obj.z(cumulation)		= obj.z(cumulation-1) + obj.vertical_velocity(cumulation-1)*dt;			
			
		end
		
        % ========================
        % integrate velocity - horizontal
        % ========================
        function [obj] = integrateVelocity(obj, cumulation, dt, direction)
            
			% Analyse desired direction (since decoupled) --> i.e. horizontal or vertical 
			if direction == "Vertical"
				velocity_error          = obj.vertical_velocity_setpoint(cumulation) - obj.vertical_velocity(cumulation -1);
				% integrate velocity
				delta_velocity          = velocity_error*obj.vertVel_p;
				obj.vertical_velocity(cumulation) = min(obj.vertical_velocity(cumulation-1)+ delta_velocity*dt, obj.max_velocity);
			elseif direction == "Horizontal"							
				% velocity error = setpoint - current velocity
				velocity_error          = obj.velocity_setpoint(cumulation) - obj.velocity(cumulation-1);
				% integrate velocity
				delta_velocity          = velocity_error*obj.vel_p;
				obj.velocity(cumulation)= min(obj.velocity(cumulation-1)+delta_velocity*dt,obj.max_velocity);
			end
		end
        		
        % ========================
        % integrate heading
        % ========================
        function [obj] = integrateHeading(obj, cumulation, dt)
            % heading error = setpoint - current heading
            heading_error           = obj.heading_setpoint-obj.heading(cumulation-1);
            % constrain heading error to the maximum yaw rate according to direction
            if heading_error >= 0
                delta_yaw           = min(obj.yaw_p*heading_error,obj.max_yawrate);
            else
                delta_yaw           = max(obj.yaw_p*heading_error,-obj.max_yawrate);
            end
            obj.heading(cumulation) = wrapTo2Pi(obj.heading(cumulation-1)+delta_yaw*dt);
        end
        
		
		% ========================
		% update target energy reserve
		% ========================
		function [obj] = updateEnergy(obj, cumulation, dt)
			% energy expenditure for moving / energy gain for successful
			% foraging
			% Assume energy expenditure scales with locomotion mode and
			% distance flown/acceleration change in modifying heading or
			% speed; 

% 			BMR = 0.877; % W, Basal metabolic rate of European starling
% 			m_bird = 0.079; % Average mass in kg of a common starling; note, they typically vary between [58, 100]g (could try and randomise between this interval)
% 			P_met_burst = 250.05*obj.m_bird^(0.8741); % Metabolic costs associated with short bursts of flight, valid for birds up to 150g in weight - sourced from Nudds and Bryant
			
% 			P_met_burst = 25; % W, average power requirement for very short
% 			burst flights (<~1 min) (for European starling) --> look for
% 			reference to help further justify this range of values (Rose
% 			Coloured Starlings, Engel S)
% 			P_met_long = 10; % W, average power requirement associated for short flights (~<10 min) European Starling

			% Check acceleration of target to determine energy cost
			% associated tith action
			if (obj.velocity(cumulation) - obj.velocity(cumulation-1))/dt >= 1.5 % m/s/s We only care about high +ve acceleration in contributing to higher power costs
				P_met = obj.P_met_burst;
			else
				P_met = obj.P_met_long;
			end

			dE = - P_met*dt; % Change in energy, ONLY WORKS WHILE dt = 1 sec increments! % P = E/t 

			% Check whether target is stationary (indicates successful
			% foraging location has been reached) -> existence energy
			if obj.velocity(cumulation-1) < 0.3 % 0.3 m/s is a fairly arbitrary but small amount
				dE = -obj.BMR*dt; % existence energy, begin decreasing at much smaller rate
				%dE = P_met*dt; % begin GAINING energy, same 'costs' as flying atm...
			end

			% distance_change = sqrt((obj.x(cumulation)-obj.x(cumulation-1))^2+(obj.y(cumulation)-obj.y(cumulation-1))^2);
			% If bird has left property stop calculating new energy
			if ~obj.in_property
				obj.energy(cumulation) = obj.energy(cumulation-1);
			else
				% change energy while still in property
				obj.energy(cumulation) = obj.energy(cumulation-1) + dE; %flying_energy*distance_change;
			end
			
			% If this flock has run out of energy it should leave property
			% now
			if obj.energy(cumulation) <= 0
				obj.in_property = false;
			end
			
			
			
			
			
		end
		       
        
        %% private methods
        % ========================
        % expand arrays
        % ========================
        function [obj] = expandArrays(obj, world, iter)
            obj.x(2:iter,:)         = ones(iter-1,1);
            obj.y(2:iter,:)         = ones(iter-1,1);
            obj.z(2:iter,:)         = ones(iter-1,1)*obj.z(1);
            obj.heading(2:iter,:)   = ones(iter-1,1);
            obj.velocity(2:iter,:)  = ones(iter-1,1);
			obj.vertical_velocity(2:iter,:)  = ones(iter-1,1);
            obj.flight(2:iter,:)    = ones(iter-1,1);
			obj.energy(2:iter,:)	= ones(iter-1,1);
            
            obj.x_setpoint(2:iter,:)= ones(iter-1,1);
            obj.y_setpoint(2:iter,:)= ones(iter-1,1);
            obj.z_setpoint(2:iter,:)= ones(iter-1,1);
            
            obj.inter_map(:,:,2:iter) = ones(world.number_of_tiles_y, world.number_of_tiles_x, iter-1);
        end
    end
end