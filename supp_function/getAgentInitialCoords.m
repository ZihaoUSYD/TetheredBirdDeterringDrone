function grid_node_coords = getAgentInitialCoords(numAgents, property_size, buffer_size, sim_num)

width = property_size(1);
height = property_size(2);

if numAgents == 1
	x = width/2;
	y = height/2;
	grid_node_coords = [x,y];
end

if numAgents == 2
	x_mid = width/2;
	y_mid = height/2;
	if width >= height
		x1 = 0.5 * x_mid;
		y1 = y_mid;
		x2 = 1.5 * x_mid;
		y2 = y_mid;
	elseif width < height
		x1 = x_mid;
		y1 = 0.5*y_mid;
		x2 = x_mid;
		y2 = 1.5*y_mid;
	else
		% It is a square..
		% Could try and implement later to account for the several unique
		% configurations that there are for 2 agents...
	end
	
	grid_node_coords = [x1 y1; x2 y2];
	
end

if numAgents > 2
	grid_node_coords = cvt_2d_sampling(numAgents,50,5000,height,width, sim_num);
end

% Amend coordinates by incorporating translation/shift from buffer region
grid_node_coords(:,1) = grid_node_coords(:,1) + buffer_size(1)/2; % X coordinate shift
grid_node_coords(:,2) = grid_node_coords(:,2) + buffer_size(2)/2; % Y coordinate shift


end