function active_agents = createActiveAgents(numActiveAgentsActual, activeAgentsCoords, maxVelocity)
% This function creates and returns an array of 'active agent' Aircraft objects
% activeAgentsCoords = createGrid(noActiveAgents,property_size(1),property_size(2));
% USING NEW CVT METHOD - note that Delaunay triangles require at least 3
% separate points to define the CVT (hence for cases of 1-2 agents will
% need to define special cases)
% activeAgentsCoords = getAgentInitialCoords(noActiveAgents,property_size(1),property_size(2), sim);

active_agents = Aircraft.empty(numActiveAgentsActual,0);
% Generate agent Aircraft structures at the desired coordinates, all
	% pointing in the same direction
	for i = 1:numActiveAgentsActual
		x_temp = activeAgentsCoords(i,1);
		y_temp = activeAgentsCoords(i,2);
		z_temp = 1.5;%1.5; % Assumes all agents must start on the ground! (i.e. before simulation takes place)
		
		% %                              typ  x     y		z		ro pc
		% yaw        ve  max_ve  min_ve	max_yawrate gain_array
		current_agent =				Aircraft(2, x_temp, y_temp, z_temp, 0, 0, deg2rad(0),    0, maxVelocity,     0,      deg2rad(90));

		% Save into agent array
		active_agents(i,1) = current_agent;
	end

	%%%% TROUBLESHOOTING CODE ---- TESTING MANUAL ENTRIES OF AGENTS -----
% 	bonus_agent = [Aircraft(2, 400, 400, z_temp, 0, 0, deg2rad(0),    0,  10,     0,      deg2rad(45))];% ...
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

	%%%%% END TROUBLESHOOTING FOR MANUAL ACTIVE AGENT ENTRIES --------
	




end