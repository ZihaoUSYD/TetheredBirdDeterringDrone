function agents_GC = createGCAgents(world,buffer_size,GC_angles)
%% This function creates and returns an array of 'ground camera agent' Aircraft objects
% Input
% world = world object
% world_buffer = buffer region outside property (for plotting purposes)
% GC_angles = 1x4 array of desired ground camera angle orientations
% (starting from lower left corner, moving C.W.)

	% Long-range camera sensors (stationary) -- one camera at each corner
	% Initialise ground camera sensors (stationary agents = Type 3)
	agents_GC = []; 

	% Extract X and Y coordinates of the world centre
	x_centre = world.world_centre(1);
	y_centre = world.world_centre(2);

	
	% Signs to find relative corner positions of a rectangle from its centre (C.W. from bottom left)
	signs_x = [-1,-1,+1,+1];
	signs_y = [-1,+1,+1,-1];


	% Create all the ground camera agents and save into array (only location
	% and pointing angles will vary between cameras)
	for i = 1:4
		temp_x = x_centre + signs_x(i)*(world.size_x-buffer_size(1))/2;
		temp_y = y_centre + signs_y(i)*(world.size_y-buffer_size(2))/2;
		temp_z = 0;

		agents_GC = [agents_GC; Aircraft(3, temp_x , temp_y, temp_z, 0, 0, deg2rad(GC_angles(i)), 0, 0, 0, 0)];
	end
end
