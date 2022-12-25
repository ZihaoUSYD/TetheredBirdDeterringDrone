function [handle] = drawRectangleColor(centre, edge_length_x, edge_length_y, boundaryColor, lineThickness)

x                   = centre(1);
y                   = centre(2);

x_bottom_l          = x-edge_length_x/2;
y_bottom_l          = y-edge_length_y/2;

x_bottom_r          = x+edge_length_x/2;
y_bottom_r          = y-edge_length_y/2;

x_top_r             = x+edge_length_x/2;
y_top_r             = y+edge_length_y/2;

x_top_l             = x-edge_length_x/2;
y_top_l             = y+edge_length_y/2;

hold on;
plot3([x_bottom_l,x_bottom_r],[y_bottom_l,y_bottom_r],[2,2],'Color', boundaryColor,'LineStyle','-','linewidth',lineThickness);
plot3([x_bottom_r,x_top_r],[y_bottom_r,y_top_r],[2,2],'Color', boundaryColor,'LineStyle','-','linewidth', lineThickness);
plot3([x_top_r,x_top_l],[y_top_r,y_top_l],[2,2],'Color', boundaryColor,'LineStyle','-','linewidth',lineThickness);
handle = plot3([x_top_l,x_bottom_l],[y_top_l,y_bottom_l],[2,2],'Color', boundaryColor,'LineStyle','-','linewidth',lineThickness);