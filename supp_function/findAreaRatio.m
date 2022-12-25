function [areaCoverageRatio,areaPropertyCoverage] = findAreaRatio(world,property_size, agents, tether_length)

world_centre = world.world_centre;

% Check validity of tether length (i.e. must be > 0)
if tether_length <= 0
	areaCoverageRatio = 0;
	areaPropertyCoverage = 0;
	return
end


%% Use approximate polygons and take the union of all the agent individual
% areas to form the total possible area coverage (accounting for overlap)

% Extract individual agent coordinates and incorporate uniform tether
% length

% First agent
agentShape = nsidedpoly(1000,'Center', [agents(1).x_tether_point, agents(1).y_tether_point], 'Radius', tether_length); % Forms a circle
propertyShape = nsidedpoly(4, 'Center', world_centre, 'SideLength', property_size(1)); % Forms a square polygon

% Now add successive agent coverage areas through union function
for i = 2:size(agents(:,1))
	tmp_agent = agents(i); % will need this later when accounting for overlapping areas
	newAgentShape = nsidedpoly(1000,'Center', [tmp_agent.x_tether_point, tmp_agent.y_tether_point], 'Radius', tether_length);
	
	% Take the union of the current agentShape amalgam with the new added
	% circle
	agentShape = union(agentShape, newAgentShape);
end

% Find the total shape of the agents which validly lies within the square
% property. I.e. take intersection of shapes
intersectedShape = intersect(agentShape, propertyShape);

% Find the area of the total agent coverage within the property
areaPropertyCoverage = area(intersectedShape);

% Find the area coverage ratio w.r.t. to property area
areaProperty = area(propertyShape);
areaCoverageRatio = areaPropertyCoverage / areaProperty;

end


