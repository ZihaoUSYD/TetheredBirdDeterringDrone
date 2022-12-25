function targets = createTargets(world,buffer_size,targetBird, numTargets)
% This function creates and returns an array of 'target' Aircraft objects
	% Initialise targets array to populate
	targets = []; 

	% Generate target parameters based off bird type

	
	
% 	BMR = 0.877; % W, Basal metabolic rate of European starling
% 	m_bird = 0.079; % Average mass in kg of a common starling; note, they typically vary between [58, 100]g (could try and randomise between this interval)
% 	P_met_burst = 250.05*m_bird^(0.8741);
% 	P_met_long = 10; % W, average power requirement associated for short flights (~<10 min) European Starling
	

	% Extract X and Y coordinates of the world centre
	x_centre = world.world_centre(1);
	y_centre = world.world_centre(2);
	
	% initialise target (make the target begin at the right end of the
	% property)
% 	target_x = x_centre + (world.size_x-buffer_size(1))/2; %  + 2*world.grid_size;
% 	target_y = y_centre - (world.size_x-buffer_size(1))/6;

	
	% Randomize target entry position (could be anywhere along the property
	% boundary, also assume that the initial heading of the flock is within a 20 degree
	% cone of the centre - somewhat arbitrary but allows for further random
	% variation)
	% Define the coordinates of the wall functions
	right_wall = x_centre + (world.size_x-buffer_size(1))/2; % x along right wall
	left_wall = x_centre - (world.size_x-buffer_size(1))/2; % x along left wall
	top_wall = y_centre + (world.size_x-buffer_size(1))/2; % y along top wall
	bottom_wall = y_centre - (world.size_x-buffer_size(1))/2; % y along bottom wall
	propertySidesFunctions = [right_wall, left_wall, top_wall, bottom_wall];
	wallOptions = [1,2,3,4]; % in order from right, top, left, bottom
	
	% Allocate start coordinates to targets
	for i = 1:numTargets
		% Randomly sample which wall to enter from
		wallChoice = randsample(wallOptions,1);
		
		% Randomly allocate remaining coordinate
		if wallChoice == 1 % right wall
			target_x_random = right_wall;
			target_y_random = bottom_wall + (top_wall-bottom_wall)*rand();
		elseif wallChoice == 2 % top wall
			target_y_random = top_wall;
			target_x_random = left_wall + (right_wall - left_wall)*rand();
		elseif wallChoice == 3 % left wall
			target_x_random = left_wall;
			target_y_random = bottom_wall + (top_wall-bottom_wall)*rand();
		elseif wallChoice == 4 % bottom wall
			target_y_random = bottom_wall;
			target_x_random = left_wall + (right_wall - left_wall)*rand();
		end
			
		% Allocate random heading buffer
		heading_buffer = deg2rad(0);
		
% 		if i == 1
% 			target_x_random = 38.6;
% 			target_y_random = 120;
% 		elseif i == 2
% 			target_x_random = 120;
% 			target_y_random = 28;
% 		else
% 			target_x_random = 120;
% 			target_y_random = 97.8;			
% 		end

% 		
		% required heading to be pointed at centre from starting position
		heading_to_centre = atan2(y_centre-target_y_random,x_centre-target_x_random);
		
		% Add noise to heading
		heading_lower_bound = heading_to_centre - heading_buffer/2;
		heading_upper_bound = heading_to_centre + heading_buffer/2;
		heading_modified = heading_lower_bound + (heading_upper_bound - heading_lower_bound)*rand();
		


		% Currently, target initial coordinates and relevant flight parameters
		% are hard-coded (for position) and assumed based off personal
		% experience (for velocity etc.)
		%                            flight		x				y			z		heading     intial_vel snse_range  max_vel	yaw_rate		target bird type/name%
		target             = Target(true, target_x_random, target_y_random,  15, heading_modified, 2	,  20,			10	,      deg2rad(20), targetBird);
		
		% Add current target to array
		targets = [targets; target];
	
	end
	
% 	target2             = Target(true, 250,			400,	20, deg2rad(-30),2,  20,			10,			deg2rad(20), targetBird);
% 	target1             = Target(true, 300, 600,  20, deg2rad(-70), 8	,  200,			10	,      deg2rad(20));

% 	targets             = [target1];%; target2];

end