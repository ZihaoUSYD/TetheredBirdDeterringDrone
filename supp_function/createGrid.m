function agentGridCoords = createGrid(numAgents,width,height)
% ===================================================================
% This function creates a grid of active agents to fill a given area
% (depends on number of desired agents, as well as property area size)
% Input:
%		numAgents		=		number of desired active agents
%		length			=		length of rectangular property
%		width			=		width of rectangular property

% Output:
%		agentGridCoords =		nx2 matrix of grid x,y coordinates representing all
%		active agents (where n = number of agents)

% Created by: Joshua Trethowan
% Created date: 8/4/2021
% ===================================================================

rectangleArea	= width*height;
pointArea		= rectangleArea/numAgents;
aspectRatio		= width/height;

% Approximate individual point areas as smaller rectangles of identical
% aspect ratio
pointDiam = sqrt(4*pointArea/pi()); % diameter between adjacent points
pointHeight = pointDiam/aspectRatio;

% Introduce buffers on height and width so that the grid is more centrally
% located in rectangle
max_agents_vert = floor(height/pointHeight);
max_agents_horz = floor(width/pointDiam);
height_buff = (height - max_agents_vert*pointHeight)/2;
width_buff  = (width - max_agents_horz*pointDiam)/2;

cols = width_buff:pointDiam:width-width_buff;
rows = height_buff:pointHeight:height - height_buff;

% Pre-create matrix
agentGridCoords = NaN(numAgents,2);
agentCounter = 1;

for i = 1:length(cols)
	for j = 1:length(rows)
		
		x_coord = cols(i);
		y_coord = rows(j);
		
		agentGridCoords(agentCounter,:) = [x_coord,y_coord];
		agentCounter = agentCounter + 1;
% 		plot(i,j,'k.','MarkerSize',14);
% 		hold on
% 		counter = counter + 1;
	end
end

% Remove any trailing NaNs
agentGridCoords(agentCounter:end,:) = [];


end