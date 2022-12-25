% ===================================================================
% This script draws the simulation environment
% void
%
% Created by: Zihao Wang
% Created date: 4 July 2017
% ===================================================================
switch world.shape
    case 'circle'
        [propH] = drawCircle(world.world_centre, world.size_other);
    case 'square'
        [propH] = drawSquare(world.world_centre, world.size_other);
    case 'rectangle'
        [propH] = drawRectangle(world.world_centre, world.size_x, world.size_y);
    otherwise
        warning('Invalid simulation environment shape argument, defaults to square.');
        [propH] = drawSquare(world.world_centre, world.world_size);
end

axis equal; axis tight; grid on;