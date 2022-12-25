function obj = bufferWaypointOLD(obj,dist_buff,yaw_buff,x_set_old,y_set_old,agent_x,agent_y,cumulation)
%% Buffer the waypoint if the next one is too close to the current one and 
% will cause unnecessary yaw oscilattions about a point (i.e. only allow
% new waypoint if it is sufficiently far or is a small angle rotation)

distance_change = euclideanDistance(x_set_old,y_set_old,agent_x,agent_y);
yaw_change = wrapTo2Pi(atan2((y_set_old-agent_y),(x_set_old-agent_x)));

if yaw_change > yaw_buff && distance_change < dist_buff %abs(new_x_setpoint - agent_x) < osc_buff || abs(new_y_setpoint - agent_y) < osc_buff
	% Keep old waypoint
	obj.x_setpoint(cumulation) = agent_x;
	obj.y_setpoint(cumulation) = agent_y;
else
	% New waypoint is sufficiently far and will not have
	% too large of a heading change
	obj.x_setpoint(cumulation) = new_x_setpoint;
	obj.y_setpoint(cumulation) = new_y_setpoint;
end

end