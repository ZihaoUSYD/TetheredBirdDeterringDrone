function figure_handle = plotTether(agent_coords,tether_coords)
x1 = agent_coords(1);
y1 = agent_coords(2);
z1 = agent_coords(3);

x2 = tether_coords(1);
y2 = tether_coords(2);
z2 = tether_coords(3);

if x1 == x2 && y1 == y2 && z1 == x2
	% Coincident/At the tether point, there is no extended cable
	return
else
	figure_handle = plot3([x1,x2],[y1,y2],[z1,z2],'k','linewidth',3);
end

end