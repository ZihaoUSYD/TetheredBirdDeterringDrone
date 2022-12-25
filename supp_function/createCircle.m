function xy_coords = createCircle(x0,y0,r)
% ===============================================================
% This generates the perimeter for a 2D circle with given coordinates and
% radius
%
% Created     : 25 March 2021
% Last updated: 25 March 2021
% By: Joshua Trethowan
% ===============================================================
theta = 0:pi()/200:2*pi();

x = x0 + r*cos(theta);
y = y0 + r*sin(theta);

xy_coords = [x;y];
end