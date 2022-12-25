% ===================================================================
% This is the main function of the UAV surveillance project
%
% agent separation rad = 130
%
% Created by: Zihao Wang
% Created date: 4 July 2017
% Last updated: 14 April 2021 (by Joshua Trethowan)
% ===================================================================
clc; clear; close all;
animation_f         = false;
videoWrite_f        = false;
addpath('supp_function');
addpath('plotting');

if videoWrite_f
    video_saver     = VideoWriter('simulation.avi');
    video_saver.FrameRate = 5;
    open(video_saver);
end

%% Performance parameters to cycle through
% No. agents
numAgentsArray = 4;
% Tether radius
tetherRadiusRange = 60;%:50:250;
% max UAV flight velocity
maxAgentSpeed = 5:2.5:20;

% Sim duration
t_sim = 150; % seconds

% Results to save
numActualAgents		= NaN(length(numAgentsArray),1);
wasSuccessful		= NaN(length(numAgentsArray),1);
timeRequired		= NaN(length(numAgentsArray),1);
tetherRadiusHor		= NaN(length(numAgentsArray),1);
agentMaxAllowSpeed	= NaN(length(numAgentsArray),1);

% Combine to form total results matrix
totalPerformance	= [numActualAgents, wasSuccessful, timeRequired];

%% Set up the simulation world (circle is not working)
% Property size (i.e. a subset of total world)
property_size = [400,400]; % [x, y] dimensions in m

% Buffer size (representing land outside property)
buffer_size = [400, 400]; % [x_buffer, y_buffer] in m, evenly split on both sides of respective dime

% Total world size
total_world_size = [property_size(1)+buffer_size(1), property_size(2) + buffer_size(2)];

legendStringArray = [];
for q = 1:length(tetherRadiusRange)
	% Reset performance matrix for the current tether value
	totalPerformance	= [numActualAgents, wasSuccessful, timeRequired];
	
for p = 1:length(numAgentsArray)
%                           shape         x    y     grid size
world               = World('rectangle',total_world_size, 20); % [1000+400,800+400]

%% Set up agents
noActiveAgents = numAgentsArray(p); % max number of active agents
maxHorizontalMovement = tetherRadiusRange(q); % maximum horizontal translation of active agents

% Create grid of agents
activeAgentsCoords = createGrid(noActiveAgents,property_size(1),property_size(2));

% Amend coordinates by incorporating translation/shift from buffer region
activeAgentsCoords(:,1) = activeAgentsCoords(:,1) + buffer_size(1)/2; % X coordinate shift
activeAgentsCoords(:,2) = activeAgentsCoords(:,2) + buffer_size(2)/2; % Y coordinate shift

% Skip cases where there are the same number of active agents (to speed
% up calcs) as before
if p > 1 % skip first case
	if length(activeAgentsCoords(:,1)) == totalPerformance(p-1,1)
		continue
	end
end

numActiveAgentsActual = length(activeAgentsCoords(:,1));
active_agents = Aircraft.empty(numActiveAgentsActual,0);

fprintf("\nCurrently on the case of %.d agents out of possible %.d max and tether = %.f m\n",numActiveAgentsActual,numAgentsArray(end),maxHorizontalMovement);

	% Generate agent Aircraft structures at the desired coordinates, all
	% pointing in the same direction
	for i = 1:numActiveAgentsActual
		x_temp = activeAgentsCoords(i,1);
		y_temp = activeAgentsCoords(i,2);
		z_temp = 0;
		% %                              typ  x     y		z		ro pc    yaw        ve  max_ve  min_ve	max_yawrate
		current_agent =				Aircraft(2, x_temp, y_temp, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45));

		% Save into agent array
		active_agents(i,1) = current_agent;

	end

	bonus_agent = [Aircraft(2, 400, 400, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45))];% ...
% 					Aircraft(2, 475, 315, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45));...
% 					Aircraft(2, 315, 315, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45))
% 					Aircraft(2, 315, 475, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45))];

	% % %                              typ  x     y    z   ro pc    yaw        ve  max_ve  min_ve max_yawrate
	% agent1              = Aircraft(2,   450,  400, 0, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45));
	% agent2              = Aircraft(2,   450,  800, 0, 0, 0, deg2rad(45),   0,  10,     0,      deg2rad(45));
	% agent3              = Aircraft(2,   950,  400, 0, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent4              = Aircraft(2,   950,  800, 0, 0, 0, deg2rad(270),  0,  10,     0,      deg2rad(45));
	% agent5              = Aircraft(2,   700,  600, 0, 0, 0, deg2rad(300),  0,  10,     0,      deg2rad(45));

	% agent6              = Aircraft(2,   500,  850, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent7              = Aircraft(2,   700,  315, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent8              = Aircraft(2,   700,  580, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent9              = Aircraft(2,   700,  850, 50, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45));
	% agent10             = Aircraft(2,   900,  315, 50, 0, 0, deg2rad(45),   0,  10,     0,      deg2rad(45));
	% agent11             = Aircraft(2,   900,  580, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent12             = Aircraft(2,   900,  850, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent13             = Aircraft(2,   1100,  315, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent14             = Aircraft(2,   1100,  580, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));
	% agent15             = Aircraft(2,   1100,  850, 50, 0, 0, deg2rad(140),  0,  10,     0,      deg2rad(45));

	% active_agents = [agent1; agent2; agent3; agent4; agent5]; % agent6; agent7; agent8; agent9; agent10; agent11; agent12; agent13; agent14; agent15];


	% Long-range camera sensors (stationary) -- always implement one on each
	% corner?
	agents_GC = [];

	% X and Y coordinates of the four rectangle corners
	x_centre = world.world_centre(1);
	y_centre = world.world_centre(2);

	GC_angles = [45,-45,225,135]; % pointing angles of cameras, moving clockwise from bottom left corner

	% Signs to find relative corner positions of a rectangle from its centre (C.W. from bottom left)
	signs_x = [-1,-1,+1,+1];
	signs_y = [-1,+1,+1,-1];


	% Create all the ground camera agents and save into array (only location
	% and pointing angles will vary between cameras)
	for i = 1:4
		temp_x = x_centre + signs_x(i)*(world.size_x-buffer_size(1))/2;
		temp_y = y_centre + signs_y(i)*(world.size_y-buffer_size(2))/2;
		temp_z = 50;

		agents_GC = [agents_GC; Aircraft(3, temp_x , temp_y, temp_z, 0, 0, deg2rad(GC_angles(i)), 0, 0, 0, 0)];
	end

	% agent_GC1           = Aircraft(3, 200, 200, 50, 0, 0, deg2rad(45), 0, 0, 0, 0);
	% agent_GC2           = Aircraft(3, 200, 1000, 50, 0, 0, deg2rad(-45), 0, 0, 0, 0);
	% agent_GC3           = Aircraft(3, 1200, 1000, 50, 0, 0, deg2rad(225), 0, 0, 0, 0);
	% agent_GC4           = Aircraft(3, 1200, 200, 50, 0, 0, deg2rad(135), 0, 0, 0, 0);

	% all agents              = [agent1;agent2;agent3];
	agents              = [active_agents; bonus_agent; agents_GC];

	% initialise agents with initial probability and environment uncertainty
	for i = 1:size(agents,1)
		agents(i)       = initProbability(agents(i), world);
		agents(i)       = initUncertainty(agents(i), world);
	end

%% Set up targets
	% initialise target (make the target begin at the right end of the
	% property)
	target_x = x_centre + (world.size_x-buffer_size(1))/2; %  + 2*world.grid_size;
	target_y = y_centre - (world.size_x-buffer_size(1))/6;

	%                            flight		x		y		z   heading       vel snse_range  max_vel		yaw_rate
	target1             = Target(true, target_x, target_y,  20, deg2rad(140), 5	,  200,			10	,      deg2rad(20));
	target2             = Target(true, 250,			400,	20, deg2rad(-30),5,  200,			10,			deg2rad(20));
% 	target1             = Target(true, 300, 600,  20, deg2rad(-70), 8	,  200,			10	,      deg2rad(20));

	targets             = [target1]%; target2];

	% initialise targets with initial interests
	for k = 1:size(targets,1)
		targets(k)      = initInterests(targets(k), world);
	end

%% plotting world
	aircraft_scale      = 30;
	multirotor_scale    = 80;
	target_scale        = 90;
	color_pallete       = ['r';'y';'k';];
% 	color_pallete_agent = ['b';'g';'c';'m';'c';'c';'y';'y';'y';'y'];
% 	color_pallete_agent = ['b';'g';'c';'m';'y';'c';'c';'c';'c';'c';'c';'c';'c';'c';'c';'y';'y';'y';'y'];
% 	plotWorld;

	color_pallete_agent1 = repmat('b',numActiveAgentsActual +1 ,1);
	color_pallete_agent2 = repmat('y',numActiveAgentsActual +1 ,1);
	color_pallete_agent3 = repmat('g',numActiveAgentsActual +1 ,1);
	
	color_pallete_GC	= repmat('c',4,1);
	color_pallete_agent_state1 = [color_pallete_agent1; color_pallete_GC];
	color_pallete_agent_state2 = [color_pallete_agent2; color_pallete_GC];
	color_pallete_agent_state3 = [color_pallete_agent3; color_pallete_GC];
	% 
	% color_pallete_agent(end-3:end) = ['y';'y';'y';'y']; % Adding ground cameras
	plotWorld;

%% Initialising simulation
	% start simulation
	t_init              = 0;
	t                   = t_init; % seconds
	dt                  = 1; % seconds
	t_total             = t_sim;
	cumulation          = 1;
	iterations          = ceil(t_total/dt)+1;
	cum_max				= 150; % simulates time cap/sim failed in deterring birds

	% expand matrices to speed up code
	world               = expandArrays(world, iterations);
	for i = 1:size(agents,1)
		agents(i)       = expandArrays(agents(i), world, iterations);
	end
	for k = 1:size(targets,1)
		targets(k)      = expandArrays(targets(k), world, iterations);
	end
	
	% Number of targets will define the number of clusters in k-means
	% clustering algorithm
	num_targets = size(targets,1);

	% skip initial condition before the loop starts
	cumulation          = cumulation+1;
	t                   = t+dt;

	% Bird begins in property
	targetsInProperty = true(1,size(targets,1));
	
	% Number of birds that have left the property so far
	deterred_bird_counter = 0;
	
	if targetsInProperty
		firstTimeAtBoundary = false;
	else
		firstTimeAtBoundary = true;
	end
%% RUN SIMULATION
	while t <= t_total && (cumulation < cum_max) && (any(targetsInProperty))
		% =============================================================================================
		% update target
		for k = 1:size(targets,1)
			% update interests
			targets(k)  = updateInterests(targets(k), agents, world, cumulation);

			% move target
			targets(k)  = moveTarget(targets(k), world, cumulation, dt);

			current_target = targets(k);
			% check if target is now out of property (ignore the first five
			% seconds)
			bound_buff = 10; % metres; include a boundary buffer along property which indicates that the bird has left the vineyard
			if t > 5 && ((targets(k).x(cumulation) < (x_centre -(world.size_x-buffer_size(1))/2)+bound_buff) || (targets(k).x(cumulation) > (x_centre + (world.size_x-buffer_size(1))/2)-bound_buff)) 
				deterred_bird_counter = deterred_bird_counter + 1;
				targets(k).in_property = false;
				targetsInProperty(k) = targets(k).in_property;
			elseif t > 5 && ((targets(k).y(cumulation) < (y_centre -(world.size_y-buffer_size(2))/2)+bound_buff) || (targets(k).y(cumulation) > (y_centre + (world.size_y-buffer_size(1))/2)-bound_buff))
				deterred_bird_counter = deterred_bird_counter + 1;
				targets(k).in_property = false;
				targetsInProperty(k) = targets(k).in_property;
			else
				% This accounts for the case that a bird has left the
				% property but has since returned while other birds are
				% still present
				targets(k).in_property = true;
				targetsInProperty(k) = targets(k).in_property;
			end
			
		end % end target loop

		% update global interests
		world           = updateGlobalInterests(world, targets, cumulation);

		% =============================================================================================
		% control agents
		for i = 1:size(agents,1)
			% update uncertainty map
			agents(i)   = updateUncertainty(agents(i), world, cumulation);

			% update probability map
			agents(i)   = updateProbability(agents(i), targets, world, cumulation);

			% reset other agent coordinate matrix (only ever concerned where
			% the other agents are in the current time step
			other_agents_coordinates = [];

			% move agent
			% experimental
			% get another agents coordinates
			if size(agents,1) > size(agents_GC,1) % i.e. there are other agents in addition to the stationary 4 ground cameras
				for j = 1:size(agents,1)
					% Skip the current agent (i.e. itself) and any ground
					% cameras
					if j==i || agents(j).type == 3
						continue
					else
						% Extract x and y coords of another agent and add them to
						% other agent matrix
						other_agents_coordinates_temp = [agents(j).x(cumulation-1),agents(j).y(cumulation-1)];
						other_agents_coordinates = [other_agents_coordinates; other_agents_coordinates_temp];
					end
				end
			end

			% Can only move active agents
			if agents(i).type == 2
				agents(i)   = moveAgent(agents(i), world, cumulation, dt, other_agents_coordinates, maxHorizontalMovement, num_targets);
				% delete max prob marker placed by previous agent
				if i == 1
					delete(agents(size(agents,1)).figure_handle_maxprob);
				else
					delete(agents(i-1).figure_handle_maxprob);
				end
			end

		end % end agent loop

		% update global uncertainty
		world           = updateGlobalUncertainty(world, agents, cumulation);

		% update global probability
		world           = updateGlobalProbability(world, agents, cumulation);

		% =============================================================================================

%%%%%%%%%%% UNCOMMENT BELOW FOR PLOTTING
		% ++++++++++++++++++++ debug figure +++++++++++++++++++	
		figure(h2);
	%     subplot(1,2,1);
		for j = 1:size(agents,1)
			delete(agents(j).figure_handle_debug);
			
			% Update colour flags for agent states
% 			if agents(j).state == "CHASE"
% 				color_pallete_agent = color_pallete_agent_state1;
% 			elseif agents(j).state == "STANDBY"
% 				color_pallete_agent = color_pallete_agent_state2;
% 			else
% 				% Either it's "GO_HOME" or it's a ground camera which will
% 				% automatically take the correct colour since it is the
% 				% same for every state
% 				color_pallete_agent = color_pallete_agent_state3;
% 			end	
			color_pallete_agent = getColorPalette(agents(j).state,"CHASE","STANDBY",color_pallete_agent_state1,color_pallete_agent_state2,color_pallete_agent_state3);
			
			switch agents(j).type
				case 1
					agents(j).figure_handle_debug = uav(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
							 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
							 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', aircraft_scale, 'color', color_pallete_agent(j),...
							 'wing', color_pallete_agent(j), 'linestyle', 'none');
				case 2
					agents(j).figure_handle_debug = quadrotor(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
							 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
							 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', multirotor_scale, 'body', color_pallete_agent(j),...
							 'boom', color_pallete_agent(j), 'prop', color_pallete_agent(j), 'linestyle', 'none');
				otherwise
					agents(j).figure_handle_debug = uav(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
							 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
							 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', aircraft_scale, 'color', color_pallete_agent(j),...
							 'wing', color_pallete_agent(j), 'linestyle', 'none');
			end
			% draw position setpoint
			delete(agents(j).figure_handle_waypoint);
			
			agents(j).figure_handle_waypoint = plot3(agents(j).x_setpoint(cumulation), agents(j).y_setpoint(cumulation), agents(j).z(cumulation), ...
					'o', 'markersize', 10, 'color', color_pallete_agent(j),'linewidth',3);
		end
		for m = 1:size(targets,1)
			delete(targets(m).figure_handle_debug);
			targets(m).figure_handle_debug = birds(targets(m).x(cumulation), targets(m).y(cumulation), targets(m).z(cumulation), ...
						'body', color_pallete(m), 'yaw', rad2deg(targets(m).heading(cumulation)), 'scale', target_scale);
			delete(targets(m).figure_handle_waypoint);
% 	        targets(m).figure_handle_waypoint = plot3(targets(m).x_setpoint(cumulation), targets(m).y_setpoint(cumulation), targets(m).z(cumulation), ...
% 	                    '^', 'markersize', 10, 'color', color_pallete(m));
		end
		delete(debug_map_handle);
		debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.prob_global(:,:,cumulation-1), 'facealpha', 1, 'linestyle', '-');
	%     debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.inter_global(:,:,cumulation-1), 'facealpha', 1, 'linestyle', '-');
% 	    debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,agents(1).cost, 'facealpha', 1, 'linestyle', '-');
		title(['Probability map - time=',num2str(t),'s'],'FontSize',30,'FontName','Arial','FontWeight','Bold');
		colorbar;

	%     subplot(1,2,2);
	%     plot(cumulation, sum(sum(world.prob_global(:,:,cumulation-1)))/world.number_of_tiles,'--or','markersize',10,'linewidth',2);
	%     hold on;
	%     title('Normalised global probability vs Time');
	%     xlabel('Time (s)'); ylabel('Global probability');
		pbaspect([1 1 1]);
		set(gcf, 'Position', [100 100 1200 800]);
		set(gca, 'FontSize', 12);

		% plot trails
		for i = 1:size(agents,1)
			color_pallete_agent = getColorPalette(agents(i).state,"CHASE","STANDBY",color_pallete_agent_state1,color_pallete_agent_state2,color_pallete_agent_state3);

			linestyle_temp = strcat(color_pallete_agent(i),'o');
			plot3(agents(i).x(cumulation), agents(i).y(cumulation), agents(i).z(cumulation), linestyle_temp,'linewidth',3); hold on;      
		end

		for m = 1:size(targets,1)
			plot3(targets(m).x(cumulation), targets(m).y(cumulation), 30, 'ro','linewidth',3); hold on;
			
		end
		% plot boundaries
		drawRectangleColor(world.world_centre, world.size_x-buffer_size(1), world.size_y-buffer_size(2));

%%% ^^^^^^ UNCOMMENT ABOVE FOR PLOTTING

	    if cumulation == 2
	        % plot maximum range of each agent as reference, only need to do once
	        max_r = tetherRadiusRange(q); % HARD-CODED ATM, needs to match value in Aircraft
	
	        for i = 1:size(agents,1)-4
	            x_0 = agents(i).x_tether_point;
	            y_0 = agents(i).y_tether_point;
	            xy_coords = createCircle(x_0,y_0,max_r);
	            plot3(xy_coords(1,:),xy_coords(2,:),30*ones(length(xy_coords(1,:))),'y--','linewidth',3); 
	            hold on;
	        end
	    end

		if videoWrite_f
			current_frame = getframe(gcf);
			writeVideo(video_saver, current_frame);
		end


		% ******************** debug figure ********************
		if mod(t,10) == 0
			fprintf("Currently at t = %.d seconds\n",t);
		end
		
		% increment time step
		t               = t + dt;
		cumulation      = cumulation + 1;
		

		
	end % end simulation loop
	
	if cumulation == cum_max
		disp('Max iterations reached');
		t = t_total;
	end

	totalPerformance(p,:)	= [numActiveAgentsActual, sum(~targetsInProperty), t];
end % end looping through variation in number of agents for current tether radius

	%% CLEAN UP DATA
	% Clean up matrix and remove NaNs
	nansInPerformance = isnan(totalPerformance);

	for i = length(totalPerformance(:,1)):-1:1
		if nansInPerformance(i,1) == 1
			totalPerformance(i,:) = [];
		end
	end


	figure(100);
	plot(totalPerformance(:,1),totalPerformance(:,3));
	hold on
	xlabel('Number of active agents in grid');
	ylabel('Time taken to repel target flock (s)');
	legendString = sprintf("tether radius = %.0f m",maxHorizontalMovement);
	legendStringArray = [legendStringArray; legendString];


end % end looping through tether radius sizes (i.e. max horizontal movements)

figure(100);
legend(legendStringArray);

if videoWrite_f
    close(video_saver);
end
%animation;