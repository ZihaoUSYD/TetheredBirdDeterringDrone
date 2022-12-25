% ===================================================================
% This is the main function of the UAV surveillance project
%
% agent separation rad = 130
%
% Created by: Zihao Wang
% Created date: 4 July 2017
% Last updated: 7 July 2021 (by Joshua Trethowan)
% ===================================================================
clc; clear; clf reset; close all;
animation_f         = false;
videoWrite_f        = false;
gainTestVideos		= false;
plotting_desired	= true; 
saveConfigs			= false;
saveAllSimData		= false;
configFilename		= 'successfulConfigs.xlsx';
simDataFilename		= 'simDataCombined.xlsx';

addpath('supp_function');
addpath('plotting');

if saveConfigs
	recycle on
	delete('successfulConfigs.xlsx'); % delete old file and prepare for new x,y, config data to be written to it
end

if saveAllSimData
	recycle on
	delete('simDataCombined.xlsx'); 
end


launchedExcel		= false;


if videoWrite_f
    video_saver     = VideoWriter('simulation.avi');
    video_saver.FrameRate = 5;
    open(video_saver);
end

if gainTestVideos
	video_saver_gains = VideoWriter('gainsVideo.avi');
	video_saver_gains.FrameRate = 5;
	open(video_saver_gains);
	
end

% Changing colors
mydarkblue = [0,0,102]/255;
mylightblue = [204, 255, 255]/255;

% Agent State Colors
myCyan			= [0, 255, 255]/255; % Ground camera agents
myBlue			= [51, 51, 255]/255; % Chase state
myYellow		= [255, 255, 102]/255; % Standby
myLightGreen	= [0, 204, 0]/255; % Go home
myDarkGreen		= [0, 153, 76]/255; % At home 

% Colors for energy
myDarkOrange = [128, 64, 0]/255;
myLightOrange = [255, 204, 153]/255;

%% PERFORMANCE PARAMETERS TO CYCLE THROUGH + SIMULATION CONDITIONS
% Performance parameters	
	% No. agents
	% 	totalActive = 10;
	numAgentsArray	= 4;%1:10;
	numTargets		= 1;
	% Tether radius
	tetherRadiusRange = 15;%0:20;
	% max UAV flight velocity
	maxAgentSpeed = 12;
	% Repeating sim to account for varying grid configurations + randomised target entry points
	num_sims = 1;%40; 
	
	
	myColorGradient = [linspace(mydarkblue(1),mylightblue(1),length(numAgentsArray))', linspace(mydarkblue(2),mylightblue(2),length(numAgentsArray))',linspace(mydarkblue(3),mylightblue(3),length(numAgentsArray))'];
	% repeat simulation at every agent number to account for varying grid
	% formations due to non-unique centroidal Voronoi tessellations
	
	% Variable parameter of interest
	target_param = maxAgentSpeed;
	
% Simulation initial parameters
	t_total             = 150;
	cum_max				= 150; % simulates time cap/sim failed in deterring birds
	dt                  = 1; % seconds
	iterations          = ceil(t_total/dt)+1;
	maxStationaryTimer  = 15; % stop simulation if there is no movement in agents/targets for this duration of time (will appear as t_max on final plot)
	t_init              = 0;

	

% Results to save
numActualAgents		= NaN(length(numAgentsArray),1);
agentMaxAllowSpeed	= NaN(length(numAgentsArray),1);
wasSuccessful		= NaN(length(tetherRadiusRange),1);
timeRequired		= NaN(length(tetherRadiusRange),1);
tetherRadiusHor		= NaN(length(tetherRadiusRange),1);
energyReserve		= NaN(length(tetherRadiusRange),1);

totalPropertyCoverageFill = NaN(length(tetherRadiusRange),1);
totalAreaRatioFill		= NaN(length(tetherRadiusRange),1);
noActiveAgentsFill		= NaN(length(tetherRadiusRange),1);
maxHorizontalMovementFill = NaN(length(tetherRadiusRange),1);
agentAreaRatioFill		= NaN(length(tetherRadiusRange),1);

% Combine to form total results matrix
% totalPerformance	= [numActualAgents, wasSuccessful, timeRequired];
totalPerformance	= [tetherRadiusHor, wasSuccessful, timeRequired, energyReserve];


% Area ratio parameters
% areaRatio			= NaN(length(numAgentsArray)*length(tetherRadiusRange),1);
% areaCoverageData	= [areaRatio, areaRatio, areaRatio]; % Will fill in with [areaRatio, numAgents, tetherLength] to represent all area info in a given deterrence event
areaCoverageData = [totalPropertyCoverageFill, totalAreaRatioFill, noActiveAgentsFill, maxHorizontalMovementFill, agentAreaRatioFill];


%% SETTING UP SIMULATION WORLD
%% Set up the simulation world (circle is not working)
% Property size (i.e. a subset of total world)
property_size = [100,100]; % [x, y] dimensions in m

% Buffer size (representing land outside property)
buffer_size = [40, 40]; % [x_buffer, y_buffer] in m, evenly split on both sides of respective dimension

% Total world size
total_world_size = [property_size(1)+buffer_size(1), property_size(2) + buffer_size(2)];

%% Create world
%                           shape         x    y   , buffer,   grid size
world               = World('rectangle',total_world_size, buffer_size, 5); % [1000+400,800+400]
x_centre = world.world_centre(1);
y_centre = world.world_centre(2);


%% Cycle through gain loops
g_p = 0.2;%0.1:0.1:1;
g_v = 0.7;%g_p;
g_a = 0.4;%g_p;

gainSimData = [];

for PPP = 1:length(g_p)
	for VVV = 1:length(g_v)
		for AAA = 1:length(g_a)
			g_p_temp = g_p(PPP);
			g_v_temp = g_v(VVV);
			g_a_temp = g_a(AAA);
			g_array = [g_p_temp; g_v_temp; g_a_temp];

%% ************ LOOP 0: Loop through simulation repetitions **********
% =============================================================================================

for sim = 1:num_sims
	
	% Randomise target entry points for the new sim / create target flock
	%% Set up targets
	bird_type = 'Starling';
	targets = createTargets(world,buffer_size,bird_type, numTargets); % Need to manually specify each target within this function (for a given world/property)
	
	% For plotting legends and titles
	legendStringArray = [];
	titleStringArray = [];
	legendStringArray2 = [];
	titleStringArray2 = [];
	legendStringArray3 = [];
	titleStringArray3 = [];
	legendStringArray4 = [];
	titleStringArray4 = [];
	legendStringArray5 = [];
	titleStringArray5 = [];
	legendStringArray6 = [];
	titleStringArray6 = [];
	legendStringArray7 = [];
	titleStringArray7 = [];
	legendStringArray8 = [];
	titleStringArray8 = [];
	
%% *********** LOOP 1: Loop through maximum agent speed *********	
for v = 1:length(maxAgentSpeed)
	maxVelocity = maxAgentSpeed(v); 
	
%% ************ LOOP 2: Loop through number of agents ***********
% =============================================================================================
for p = 1:length(numAgentsArray)
	totalPerformance	= [tetherRadiusHor, wasSuccessful, timeRequired, energyReserve]; % reset performance metric for new number of agents
	areaCoverageData = [totalPropertyCoverageFill, totalAreaRatioFill, noActiveAgentsFill, maxHorizontalMovementFill, agentAreaRatioFill];

	
	%% Set up agents
	noActiveAgents = numAgentsArray(p); % max number of active agents in current simulation
	
	% Use CVT sampling to extract grid coordinates of agents (indicating
	% the current simulation number to track non-uniqueness)
	activeAgentsCoords = getAgentInitialCoords(noActiveAgents,property_size, buffer_size, sim);
		
	% Create grid of agents
	numActiveAgentsActual = length(activeAgentsCoords(:,1));
	active_agents = createActiveAgents(numActiveAgentsActual, activeAgentsCoords, maxVelocity, g_array);

	% Create Ground Camera Sensors (i.e. stationary agents)
	GC_angles = [45,-45,225,135]; % pointing angles of cameras, moving clockwise from bottom left corner (i.e. all pointing inward)
	GC_angles_matCam = deg2rad([45;135;225;315]); % plotCamera defines its orientation angles differently, use these for plotting camera and the above angles for actual sensing data
	agents_GC = createGCAgents(world,buffer_size,GC_angles, g_array);

	% Concatenate all agents
	agents = [active_agents; agents_GC];
	
	%% Initialise agents and targets with initial probability, environment uncertainty and interests
	for i = 1:size(agents,1)
		agents(i)       = initProbability(agents(i), world);
		agents(i)       = initUncertainty(agents(i), world);
	end
	
	for k = 1:size(targets,1)
		targets(k)      = initInterests(targets(k), world);
	end
	
%% ************ LOOP 3: Loop through tether radius range values *******
% =============================================================================================
for q = 1:length(tetherRadiusRange)
	%% Set up variable tether radius
	% Maximum allowed horizontal translation of active agents
	maxHorizontalMovement = tetherRadiusRange(q);

	% Troubleshooting code for seeing which case we are on in the simulation
	fprintf("\nSim:%d - Max v: %.f - currently on the case of %.d agents out of possible %.d max and tether = %.f m\n",sim,maxVelocity,numActiveAgentsActual,numAgentsArray(end),maxHorizontalMovement);

%% Plot world if flag has been checked
if plotting_desired == true
	aircraft_scale      = 10;
	multirotor_scale    = 30;
	target_scale        = 30;
	color_pallete       = ['r';'y';'k';];

	color_pallete_agent1 = repmat(myBlue,numActiveAgentsActual ,1);			% Chase
	color_pallete_agent2 = repmat(myYellow,numActiveAgentsActual ,1);		% Standby
	color_pallete_agent3 = repmat(myLightGreen,numActiveAgentsActual ,1);	% Go Home
	color_pallete_agent4 = repmat(myDarkGreen,numActiveAgentsActual, 1);	% Home
	
	color_pallete_GC	= repmat(myCyan,4,1);
	color_pallete_agent_state1 = [color_pallete_agent1; color_pallete_GC];
	color_pallete_agent_state2 = [color_pallete_agent2; color_pallete_GC];
	color_pallete_agent_state3 = [color_pallete_agent3; color_pallete_GC];
	color_pallete_agent_state4 = [color_pallete_agent4; color_pallete_GC];
		
	% Plot world figure
	% Choose view aspect
	viewAspectDesired = [-1,-1,0.55]; % This is close to isometric 
	plotWorld;
	
	% Draw ground cameras
% 	figure(h2);
% 	for j = 1:size(GC_angles_matCam,1)
% 		theta1 = pi()/2; % rotates camera about x axis onto xy plane
% 		theta2 = GC_angles_matCam(j);
% 		
% 		Rx = [1     0     0;
% 			0     cos(theta1)    -sin(theta1);
% 			0     sin(theta1)     cos(theta1)];
% 		
% 		Rz = [cos(theta2)     -sin(theta2)     0;
% 			sin(theta2)     cos(theta2)    0;
% 			0     0     1];
% 		
% 		R = Rx*Rz;
% 		x = agents_GC(j).x;
% 		y = agents_GC(j).y;
% 		z = agents_GC(j).z;
% 		
% 		plotCamera('Location', [x,y,z],'Orientation',R, 'Color', myCyan,'Size', 5);
% 		
% 		hold on
% 	end
end

%% Initialising simulation
	% start simulation / reset conditions
	t                   = t_init; % seconds
	cumulation          = 1;
	
	% skip initial condition before the loop starts
	cumulation          = cumulation+1;
	t                   = t+dt;
	
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

	% Bird begins in property
	targetsInProperty = true(1,size(targets,1));
	
	% Number of birds that have left the property so far
	deterred_bird_counter = 0;
	
	if targetsInProperty
		firstTimeAtBoundary = false;
	else
		firstTimeAtBoundary = true;
	end
	
	% Assume target start with sufficient energy reserve
	sufficientTargetEnergy = true;
	
	% Stationary target timer
	allTargetsStationaryTimer = 0;
	allAgentsOnStandby = false;
	numAgentsOnStandby = 0;
	
%% RUN SIMULATION
	while t < t_total && (cumulation < cum_max) && (any(targetsInProperty)) && sufficientTargetEnergy
		% =============================================================================================
		% update target
		for k = 1:size(targets,1)
			% update interests
			targets(k)  = updateInterests(targets(k), agents, world, cumulation);

			% move target
			targets(k)  = moveTarget(targets(k), world, cumulation, dt);
			
			% update energy reserves of target
			targets(k) = updateEnergy(targets(k), cumulation, dt);

			current_target = targets(k);
			bound_buff = 0; % metres; include a boundary buffer along property which indicates that the bird has left the vineyard	
			
			% Check if current target is now out of property (ignore the first five
			% seconds)
			% Check x coordinate
			if t > 5 && ((targets(k).x(cumulation) < (x_centre -(world.size_x-buffer_size(1))/2)+bound_buff) || (targets(k).x(cumulation) > (x_centre + (world.size_x-buffer_size(1))/2)-bound_buff)) 
				deterred_bird_counter = deterred_bird_counter + 1;
				targets(k).in_property = false;
				targetsInProperty(k) = targets(k).in_property;
			% Check y coordinate
			elseif t > 5 && ((targets(k).y(cumulation) < (y_centre -(world.size_y-buffer_size(2))/2)+bound_buff) || (targets(k).y(cumulation) > (y_centre + (world.size_y-buffer_size(1))/2)-bound_buff))
				deterred_bird_counter = deterred_bird_counter + 1;
				targets(k).in_property = false;
				targetsInProperty(k) = targets(k).in_property;
			else
				% This accounts for the case that a bird has left the
				% property but has since returned while other birds are
				% still present OR the bird just hasn't left
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

%%% Plot debug figure if flag has been checked
if plotting_desired == true
%% 		++++++++++++++++++++ debug figure +++++++++++++++++++	
		figure(h2);
	%     subplot(1,2,1);
		for j = 1:size(agents,1)
			delete(agents(j).figure_handle_debug);
			
			color_pallete_agent = getColorPalette(agents(j).state,"CHASE","STANDBY","HOME",color_pallete_agent_state1,color_pallete_agent_state2,color_pallete_agent_state3, color_pallete_agent_state4);
			
			switch agents(j).type
				case 1
					agents(j).figure_handle_debug = uav(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
							 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
							 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', aircraft_scale, 'color', color_pallete_agent(j,:),...
							 'wing', color_pallete_agent(j,:), 'linestyle', 'none');
				case 2
					agents(j).figure_handle_debug = quadrotor(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
							 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
							 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', multirotor_scale, 'body', color_pallete_agent(j,:),...
							 'boom', color_pallete_agent(j,:), 'prop', color_pallete_agent(j,:), 'linestyle', 'none');

					% draw position setpoint
					delete(agents(j).figure_handle_waypoint);
% 					agents(j).figure_handle_waypoint = plot3(agents(j).x_setpoint(cumulation), agents(j).y_setpoint(cumulation), agents(j).z_tether_point, ...
% 							'o', 'markersize', 10, 'color', color_pallete_agent(j,:),'linewidth',3);
					
					% Create circle for waypoint
					tmp_waypoint_xy = createCircle(agents(j).x_setpoint(cumulation), agents(j).y_setpoint(cumulation),1);
					tmp_waypoint_x = tmp_waypoint_xy(1,:);
					tmp_waypoint_y = tmp_waypoint_xy(2,:);
					tmp_waypoint_z = agents(j).z_tether_point*ones(1,length(tmp_waypoint_x));
					agents(j).figure_handle_waypoint = plot3(tmp_waypoint_x, tmp_waypoint_y, tmp_waypoint_z, 'color', color_pallete_agent(j,:),'Marker','.','MarkerSize',10);
					
					% Make sure to match view aspect of waypoints to that
					% of the world
% 					view(viewAspectDesired);
% 						
					% draw agent tether
					delete(agents(j).figure_handle_tether);
					temp_agent_coords = [agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation)];
					temp_agent_tether_coords = [agents(j).x_tether_point, agents(j).y_tether_point, agents(j).z_tether_point];
					agents(j).figure_handle_tether = plotTether(temp_agent_coords, temp_agent_tether_coords);
					
					
				%otherwise
					% Do nothing
					
					%agents(j).figure_handle_debug = uav(agents(j).x(cumulation), agents(j).y(cumulation), agents(j).z(cumulation),...
					%		 'roll', rad2deg(agents(j).roll(cumulation)), 'pitch', rad2deg(agents(j).pitch(cumulation)),...
					%		 'yaw', rad2deg(agents(j).yaw(cumulation)), 'scale', aircraft_scale, 'color', color_pallete_agent(j,:),...
					%		 'wing', color_pallete_agent(j,:), 'linestyle', 'none');
			end
		end
		for m = 1:size(targets,1)
			delete(targets(m).figure_handle_debug);
			targets(m).figure_handle_debug = birds(targets(m).x(cumulation), targets(m).y(cumulation), targets(m).z(cumulation), ...
						'body', color_pallete(m), 'yaw', rad2deg(targets(m).heading(cumulation)), 'scale', target_scale);
			delete(targets(m).figure_handle_waypoint);
	        targets(m).figure_handle_waypoint = plot3(targets(m).x_setpoint(cumulation), targets(m).y_setpoint(cumulation), targets(m).z(cumulation), ...
	                    '^', 'markersize', 10, 'color', color_pallete(m));
		end
		delete(debug_map_handle);
		
		
		debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.prob_global(:,:,cumulation-1), 'facealpha', 1, 'linestyle', '-', 'EdgeColor',[191 191 191]/255);
% 	    debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.inter_global(:,:,cumulation-1), 'facealpha', 1, 'linestyle', '-');
% 	    debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,agents(1).cost, 'facealpha', 1, 'linestyle', '-');
		title(['Probability map - time=',num2str(t),'s'],'FontSize',30,'FontName','Arial','FontWeight','Bold');
		colorbar;
	
		% Adjust colour of text depending on whether energy reserve is increasing or
		% decreasing
% 		myDarkGreen = [0, 153, 51]/255;
% 		energyColorPalette = [myDarkGreen,'r'];
		if targets(1).energy(cumulation) >= targets(1).energy(cumulation - 1 )
			energyColor = myDarkGreen;
		else
			energyColor = 'r';
		end
		
		delete(an);
		% add an annotation showing current target energy reserve
		an = annotation('textbox',[0.4 .85 0 0],'String',['Target energy reserve: ',num2str(round(targets(1).energy(cumulation),3,'significant')),' J'],'Color',energyColor,'FontSize',12,'BackgroundColor','w','FitBoxToText','on');
		
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
			color_pallete_agent = getColorPalette(agents(i).state,"CHASE","STANDBY","HOME",color_pallete_agent_state1,color_pallete_agent_state2,color_pallete_agent_state3, color_pallete_agent_state4);

			lineColor = color_pallete_agent(i,:);
% 			plot3(agents(i).x(cumulation), agents(i).y(cumulation), agents(i).z(cumulation), 'Color', lineColor,'LineStyle','none','Marker','.','MarkerSize',25); hold on;      
			plot3(agents(i).x(cumulation-1:cumulation), agents(i).y(cumulation-1:cumulation), agents(i).z(cumulation-1:cumulation), 'Color', lineColor,'LineWidth',2.5); hold on;      
		end

		for m = 1:size(targets,1)
% 			plot3(targets(m).x(cumulation), targets(m).y(cumulation), targets(m).z(cumulation), 'r.','MarkerSize',25); hold on;
			plot3(targets(m).x(cumulation-1:cumulation), targets(m).y(cumulation-1:cumulation), targets(m).z(cumulation-1:cumulation), 'r', 'LineWidth',2.5); hold on;

			
		end
		% plot boundaries
		drawRectangleColor(world.world_centre, world.size_x-buffer_size(1), world.size_y-buffer_size(2), 'k');


%% for h2 'Debug figure' figure handle as well

	    if cumulation == 2
	        % plot maximum range of each agent as reference, only need to do once
	        max_r = tetherRadiusRange(q); 
	
	        for i = 1:size(agents,1)-4
	            x_0 = agents(i).x_tether_point;
	            y_0 = agents(i).y_tether_point;
	            xy_coords = createCircle(x_0,y_0,max_r);
	            plot3(xy_coords(1,:),xy_coords(2,:),1.5*ones(length(xy_coords(1,:))),'k.','linewidth',3); 
	            hold on;
	        end
	    end
end
%%% ^^^^^^ End debug figure plotting
		if videoWrite_f
			current_frame = getframe(gcf);
			writeVideo(video_saver, current_frame);
		end
		
		if gainTestVideos
			current_frame2 = getframe(gcf);
			writeVideo(video_saver_gains, current_frame2);
		end


		% ******************** debug figure ********************
		if mod(t,10) == 0
			fprintf("Currently at t = %.d seconds\n",t);
		end
		
		% increment time step
		t               = t + dt;
		cumulation      = cumulation + 1;
		
		if cumulation == cum_max
			disp('Max iterations reached');
			t = t_total;
		end
		
		% Check velocity of all targets --> if they are all stationary add
		% to stationary target ticker/timer, also energy reserve
		target_vel_array_temp = NaN(size(targets,1),1);
		target_energy_array_temp = NaN(size(targets,1),1);
		
		for i = 1:size(targets,1) 
			target_v = targets(i).velocity(cumulation-1);
			target_e = targets(i).energy(cumulation-1);
			target_vel_array_temp(i) = target_v;
			target_energy_array_temp(i) = target_e;
		end
		
		if max(target_vel_array_temp) < 0.3 % 0.5 m/s is fairly aribtrary but close enough to effective zero
			allTargetsStationaryTimer = allTargetsStationaryTimer + dt; % all the targets are essentially stationary
		else
			allTargetsStationaryTimer = 0; % reset timer because there is at least one agent still moving or being chased
		end
		
		% Check states of all agents
		for i = 1:size(agents,1)
			if agents(i).state == "STANDBY" || (agents(i).velocity(cumulation-1) < 0.5 && agents(i).type == 2) 
				numAgentsOnStandby = numAgentsOnStandby + 1;
			end
		end
		
		if numAgentsOnStandby == noActiveAgents
			allAgentsOnStandby = true;
		else
			allAgentsOnStandby = false;
		end
		
		numAgentsOnStandby = 0; % reset counter for next time step
		
		% All targets are stationary have been stationary for the given
		% amount of buffer time and agents are all on standby mode,
		% indicating there will be no change; end simulation now
		if allAgentsOnStandby && (allTargetsStationaryTimer >= maxStationaryTimer)
			t = t_total;
		end
		
		% Check energy reserve of target.
		if target_e <= 0
			sufficientTargetEnergy = false; % target no longer has sufficient energy reserve for flying without fleeing; mission success
		end
		
	end % end current time simulation loop
	
% 	totalPerformance(p,:)	= [numActiveAgentsActual, sum(~targetsInProperty), t];
	missionWasSuccess		= all(~targetsInProperty) || ~sufficientTargetEnergy; % Two ways in which the mission can be successful: birds are deterred or they've used their threshold total
	totalPerformance(q,:)	= [maxHorizontalMovement, missionWasSuccess , t, targets(1).energy(cumulation-1)];
	
	% Area coverage
	[totalAreaRatio, totalPropertyCoverage] = findAreaRatio(world, property_size, agents(1:end-4,:), maxHorizontalMovement);
	agentAreaRatio = noActiveAgents / totalPropertyCoverage; % Number of agents per "coverage" area (which is dependent on tether size)
	areaCoverageData(q,:) = [totalPropertyCoverage, totalAreaRatio, noActiveAgents, maxHorizontalMovement, agentAreaRatio]; % Append the Area coverage ratio, numAgents, tetherRange
	
	
	% If saving all sim data, write to file the performance measures etc
	if saveAllSimData
		simHeadings = {'Sim no.','No. Agent/s','No. Target/s','Target x entry (m)','Target y entry (m)','Tether radius (m)','Agent max speed (m/s)','Mission success?','Time to dispel (s)','Energy expended by bird in flock (J)'};
		% This will launch the excel regardless
		xlswrite(simDataFilename,simHeadings);
		
% 		if missionWasSuccess
% 			successMessage = 'y';
% 		else
% 			successMessage = 'n';
% 		end
		
		energyExpendedBird = targets(1).energy(1) - targets.energy(cumulation-1);
		x_entry = targets(1).x(1);
		y_entry = targets(1).y(1);
		% Append data under headings / under previous sim data
		simDataToAppend = {sim, numAgentsArray(p), length(targets), x_entry, y_entry, maxHorizontalMovement, maxVelocity, missionWasSuccess, totalPerformance(q,3), energyExpendedBird};
		xlsappend(simDataFilename, simDataToAppend);
	end

	
	
end % end looping through variation in tether radius


	%% CLEAN UP DATA
	% Clean up matrix and remove NaNs
	nansInPerformance = isnan(totalPerformance);

	for i = length(totalPerformance(:,1)):-1:1
		if nansInPerformance(i,1) == 1
			totalPerformance(i,:) = [];
		end
	end
	
	%% EXTRACT AND SAVE ONLY SUCCESSFUL CASES/DATA
	% For successful cases, save data to plot tether radius vs. agent/area
	% ratio
	% Extract only 'successful' sim data for current agent number
	successfulSimData = []; % Parameters such as tether radius etc.
	successfulTetherConfigs = []; % x,y, coordinates of successful configuration as well as tether radius
	
	for indx = 1:length(totalPerformance(:,1))
		if totalPerformance(indx,2) == true % corrresponds to missionWasSuccess boolean
			successfulTethersTemp = totalPerformance(indx,1);
			successfulSimData = [successfulSimData; successfulTethersTemp, areaCoverageData(indx,2), areaCoverageData(indx,3)]; % [tether radius, agent area ratio, number of agents]
			successfulTetherConfigs = [successfulTetherConfigs; successfulTethersTemp];
		end
	end
	
% 	% Plot property area coverage vs. tether radius (FOR SUCCESSFUL SIMS - if there are any)
% 	if ~isempty(successfulSimData)
% 		figure(800+10*sim + v);
% 			plot(successfulSimData(:,1),successfulSimData(:,2)*100,'Color',myColorGradient(p,:),'LineStyle','-','LineWidth',2); % Note: Multiplied by 100 to convert to %
% 			hold on
% 			xlabel('Tether radius (m)');
% 			ylabel('Agent coverage area of property (%)')
% 			legendString8 = sprintf("%d agent/s",successfulSimData(1,3));
% 			legendStringArray8 = [legendStringArray8; legendString8];
% 			titleString8 = sprintf("SUCCESSFUL %d m x %d m grid: Simulation %d, Max V = %.f m/s", property_size(1),property_size(2),sim, maxVelocity);
% 	end
	
	
	% SAVE Agent grid configuration data of successful/unsuccessful grids
	% Note: sheet no. = agent no.
	%		row corresponds to new simulation data
	%		[x,y] coords ==> [A, B] Excel columns, tether radius ==> [C]
	%		Excel column
	% Only save data if there is data to save!
	
	if ~isempty(successfulTetherConfigs) && saveConfigs == true
		if launchedExcel == true
			currentSheetNameArray = sheetnames('successfulConfigs.xlsx');
			desiredSheetName = strcat("Sheet",num2str(numAgentsArray(p)));
			sheetAlreadyExists = any(contains(currentSheetNameArray,desiredSheetName)); % Check whether sheet already exists in Excel
			
			if sheetAlreadyExists
				% Check if the sheet has been populated already
				sheetIsEmpty = isempty(xlsread('successfulConfigs.xlsx',desiredSheetName));
			else
				% If sheet doesn't already exist, sheet must also be empty
				sheetIsEmpty = true;
			end
			
			
		end

		if launchedExcel == false 
			% Need to write to excel file the first time (include headings
			% for clarity later when analysing) 
			% HEADINGS
			sheetHeading = {"Num agents on grid:", numAgentsArray(p)}; %d agent/s",numAgentsArray(p));
			xlswrite('successfulConfigs.xlsx', sheetHeading, numAgentsArray(p),'A1');
			xlswrite('successfulConfigs.xlsx', ["x","y","tether length/s (m) - from column C onward"], numAgentsArray(p),'A2');
			% FIRST ROW OF DATA (i.e. row #2 of Excel)
			xlswrite('successfulConfigs.xlsx', [activeAgentsCoords, repmat(successfulTetherConfigs',length(activeAgentsCoords(:,1)),1)], numAgentsArray(p),'A3');
			launchedExcel = true;
			
		elseif ~sheetAlreadyExists || sheetIsEmpty
			% Sheet does not already exist, print headings
			sheetHeading = {"Num agents on grid:", numAgentsArray(p)}; %d agent/s",numAgentsArray(p));
			xlswrite('successfulConfigs.xlsx', sheetHeading, numAgentsArray(p),'A1');
			xlswrite('successfulConfigs.xlsx', ["x","y","tether length/s (m) - from column C onward"], numAgentsArray(p),'A2');
			% FIRST ROW OF DATA (i.e. row #2 of Excel)
			xlswrite('successfulConfigs.xlsx', [activeAgentsCoords, repmat(successfulTetherConfigs',length(activeAgentsCoords(:,1)),1)], numAgentsArray(p),'A3');
			
		else
			% Now we can just append data to open Excel on an existing
			% sheet
			xlsappend('successfulConfigs.xlsx', [activeAgentsCoords, repmat(successfulTetherConfigs',length(activeAgentsCoords(:,1)),1)], numAgentsArray(p));
		end
	end % end saving successful grid data into Excel for current parameters
	
	
	%% PLOTTING PART 1
% 	% Plot Time vs. Tether radius (for every given agent no.)
% 	figure(100 + 10*sim + v);
% 	plot(totalPerformance(:,1),totalPerformance(:,3),'Color',myColorGradient(p,:),'LineWidth',2);
% 	hold on
% 	xlabel('Tether radius (m)');
% 	ylabel('Time taken to repel target flock (s)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray = [legendStringArray; legendString];
% 	titleString = sprintf("%d m x %d m grid: Simulation %d, Max V = %.f m/s", property_size(1),property_size(2),sim, maxVelocity);
% 	titleStringArray = [titleStringArray; titleString];
% 	
% 	plot(totalPerformance(:,1),totalPerformance(:,3));
% 	hold on
% 	xlabel('Number of active agents in grid');
% 	ylabel('Time taken to repel target flock (s)');
% % 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString = sprintf("tether = %.0f m",maxHorizontalMovement);

	% 	figure(200 + sim);
% 	plot(totalPerformance(:,3),totalPerformance(:,4),'Color',myColorGradient(p,:),'LineWidth',2);
% 	% 	plot(totalPerformance(:,3),totalPerformance(:,4),'b','LineWidth',2);
% 	hold on
% 	xlabel('Time (s)');
% 	ylabel('Energy reserve (credits/units)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString2 = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray2 = [legendStringArray2; legendString2];
% 	titleString2 = sprintf("%d m x %d m grid: Simulation %d", property_size(1),property_size(2),sim);
% 	titleStringArray2 = [titleStringArray2; titleString2];

% 	figure(300 + sim);
% 	plot(totalPerformance(:,1),totalPerformance(:,3),'Color',myColorGradient(p,:),'LineStyle','-','LineWidth',2);
% % 	plot(totalPerformance(:,1),totalPerformance(:,3),'b','LineWidth',2);
% 	hold on
% 	xlabel('Tether radius (m)');
% 	ylabel('Time taken to repel target flock (s)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString3 = sprintf("%d agent/s",numAgentsArray(p));
% 
% 	yyaxis right
% 	plot(totalPerformance(:,1),totalPerformance(:,4),'Color',myColorEnergyGradient(p,:),'LineStyle','-','LineWidth',2);
% 	ylabel('Energy reserve (credits/units)');
% 	legendString3Energy = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray3 = [legendStringArray3; legendString3; legendString3Energy];
% 	titleString3 = sprintf("%d m x %d m grid: Simulation %d", property_size(1),property_size(2),sim);
% 	titleStringArray3 = [titleStringArray3; titleString3];

	% Plot Energy expenditure reserve vs. Time (for every given agent no.)
% 	figure(400+10*sim + v);
% 	plot(0:cumulation-2,targets(1).energy(1:cumulation-1),'Color',myColorGradient(p,:),'LineWidth',2);
% 	xlabel('Time (s)');
% 	ylabel('Energy reserve (units)');
% 	legendString4 = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray4 = [legendStringArray4; legendString4];
% 	titleString4 = sprintf("%d m x %d m grid: Simulation %d, Max V = %.f m/s", property_size(1),property_size(2),sim, maxVelocity);
% 	hold on

	% Plot Energy expenditure reserve vs. Tether radius (for every given agent no.)
% 	figure(500+10*sim + v);
% 	plot(totalPerformance(:,1),totalPerformance(:,4),'Color',myColorGradient(p,:),'LineStyle','-','LineWidth',2);
% % 	plot(totalPerformance(:,1),totalPerformance(:,3),'b','LineWidth',2);
% 	hold on
% 	xlabel('Tether radius (m)');
% 	ylabel('Bird flock energy reserve (units)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString5 = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray5 = [legendStringArray5; legendString5];
% 	titleString5 = sprintf("%d m x %d m grid: Simulation %d, Max V = %.f m/s", property_size(1),property_size(2),sim, maxVelocity);
	
	% Plot Time  vs. Agent coverage area RATIO (for every given agent no.)
% 	figure(600+sim);
% 	plot(areaCoverageData(:,5),totalPerformance(:,3),'Color',myColorGradient(p,:),'LineStyle','-','LineWidth',2);
% % 	plot(totalPerformance(:,1),totalPerformance(:,3),'b','LineWidth',2);
% 	hold on
% 	xlabel('Agent density - no. agents per property coverage area (agents/m^2)');
% 	ylabel('Time (s)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString6 = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray6 = [legendStringArray6; legendString6];
% 	titleString6 = sprintf("%d m x %d m grid: Simulation %d", property_size(1),property_size(2),sim);
	
	
	% Plot  Time vs. Agent area coverage(for every given agent no.) -->
	% look for trend based on tether radius??
% 	figure(700+10*sim + v);
% 	plot(areaCoverageData(:,2)*100,totalPerformance(:,3),'Color',myColorGradient(p,:),'LineStyle','-','LineWidth',2); % Note: Multiplied by 100 to convert to %
% % 	plot(totalPerformance(:,1),totalPerformance(:,3),'b','LineWidth',2);
% 	hold on
% 	xlabel('Agent coverage area of property (%)');
% 	ylabel('Time (s)');
% 	% 	area_density = (pi()*maxHorizontalMovement^2)/(property_size(1)*property_size(2));
% 	legendString7 = sprintf("%d agent/s",numAgentsArray(p));
% 	legendStringArray7 = [legendStringArray7; legendString7];
% 	titleString7 = sprintf("%d m x %d m grid: Simulation %d, Max V = %.f m/s", property_size(1),property_size(2),sim, maxVelocity);
% 	
	
	

	
	
	
	

end % end looping through variation in agent number

%% PLOTTING PART 2 (Clean-up)
% Plot results of current simulation
% figure(100 + 10*sim + v);
% legend(legendStringArray, 'Location', 'NorthEastOutside');
% title(titleString);
% grid on;
% grid minor;

% 
% figure(200 + sim);
% legend(legendStringArray2, 'Location', 'NorthEastOutside');
% title(titleString2);
% grid on;
% grid minor;

% figure(300 + sim);
% legend(legendStringArray3, 'Location', 'NorthEastOutside');
% title(titleString3);
% grid on;
% grid minor;


% Plot energy reserve variation for target flocks
% figure(400 + 10*sim + v);
% title(sprintf('Energy reserve of target bird flock: Sim %d',sim)); 
% legend(legendStringArray4, 'Location', 'NorthEastOutside');
% grid on;
% grid minor;
% 
% figure(500 +10*sim + v);
% legend(legendStringArray5, 'Location', 'NorthEastOutside');
% title(titleString5);
% grid on;
% grid minor;

% figure(600 + sim);
% legend(legendStringArray6, 'Location', 'NorthEastOutside');
% title(titleString6);
% grid on;
% grid minor;

% figure(700 + 10*sim + v);
% legend(legendStringArray7, 'Location', 'NorthEastOutside');
% title(titleString7);
% grid on;
% grid minor;
% 
% if ~isempty(legendStringArray8)
% 	figure(800 + 10*sim + v);
% 	legend(legendStringArray8, 'Location', 'NorthEastOutside');
% 	title(titleString8);
% 	grid on;
% 	grid minor;
% end

end % end variation in agent speed

end % end simulation repeats

% Print out maximum and minimum heights across the agents
max_height = max(agents(1).z);
min_height = min(agents(1).z);

for i = 2:length(agents)
	if max(agents(i).z) > max_height
		max_height = max(agents(i).z);
	end
	
	if min(agents(i).z) < min_height
		min_height = min(agents(i).z);
	end
end

fprintf('\nMax height: %.2fm, Min height: %.2fm\n Gains: Gp = %.1f, Gv = %.1f, Ga = %.1f\n',max_height, min_height, g_p(PPP),g_v(VVV),g_a(AAA));

% Save results
gainSimData = [gainSimData; g_p(PPP),g_v(VVV),g_a(AAA), max_height, min_height];

		end
	end
end


if videoWrite_f
    close(video_saver);
end

if gainTestVideos
	close(video_saver_gains);
end
%animation;