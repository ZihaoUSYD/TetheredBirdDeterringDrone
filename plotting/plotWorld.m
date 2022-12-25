% %% plot the global uncertainty
% h1                  = figure('Name', 'World', 'NumberTitle', 'off');
% % subplot(2,3,1);
% drawBoundary;
% hold on;
% % % draw aircraft and map
% for i = 1:size(agents,1)
%     switch agents(i).type
%         case 1
%             agents(i).figure_handle_uncertainty = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%         case 2
%             agents(i).figure_handle_uncertainty = quadrotor(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', multirotor_scale, 'body', color_pallete(i),...
%                      'boom', color_pallete(i), 'prop', color_pallete(i), 'linestyle', 'none');
%         otherwise
%             agents(i).figure_handle_uncertainty = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%     end
% end
% % draw targets
% for i = 1:size(targets,1)
%     targets(i).figure_handle_uncertainty = birds(targets(i).x, targets(i).y, targets(i).z, ...
%                 'body', color_pallete(i), 'yaw', rad2deg(targets(i).heading), 'scale', target_scale);
% end
% uncertainty_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.uncert_global, 'facealpha', 1, 'linestyle', 'none');
% oldcolormap = colormap('gray');
% colormap( flipud(oldcolormap) );
% caxis([0 1]);
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;
% axis([-100 world.size_x+100 -100 world.size_y+100]);
% 
% % scatter plot of envrionment uncertainty
% subplot(2,3,4);
% plot(0, sum(sum(world.uncert_global))/world.number_of_tiles,'--or','markersize',10,'linewidth',2);
% grid on;
% title('Evolution of global envrioment uncertainty');
% hold on;
% 
% %% plot the global probability
% subplot(2,3,2);
% drawBoundary;
% hold on;
% % draw aircraft and map
% for i = 1:size(agents,1)
%     switch agents(i).type
%         case 1
%             agents(i).figure_handle_probability = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%         case 2
%             agents(i).figure_handle_probability = quadrotor(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', multirotor_scale, 'body', color_pallete(i),...
%                      'boom', color_pallete(i), 'prop', color_pallete(i), 'linestyle', 'none');
%         otherwise
%             agents(i).figure_handle_probability = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%     end
% end
% % draw targets
% for i = 1:size(targets,1)
%     targets(i).figure_handle_probability = birds(targets(i).x, targets(i).y, targets(i).z, ...
%                 'body', color_pallete(i), 'yaw', rad2deg(targets(i).heading), 'scale', target_scale);
% end
% probability_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.prob_global, 'facealpha', 1, 'linestyle', 'none');
% oldcolormap = colormap('gray');
% colormap( flipud(oldcolormap) );
% caxis([0 1]);
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;
% axis([-100 world.size_x+100 -100 world.size_y+100]);
% 
% % scatter plot of target probability
% subplot(2,3,5)
% plot(0, sum(sum(world.prob_global))/world.number_of_tiles,'--ob','markersize',10,'linewidth',2);
% grid on;
% title('Evolution of global target probability');
% hold on;
% 
% %% Plot target and interests
% subplot(2,3,3);
% drawBoundary;
% hold on;
% % draw aircraft and map
% for i = 1:size(agents,1)
%     switch agents(i).type
%         case 1
%             agents(i).figure_handle_interests = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%         case 2
%             agents(i).figure_handle_interests = quadrotor(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', multirotor_scale, 'body', color_pallete(i),...
%                      'boom', color_pallete(i), 'prop', color_pallete(i), 'linestyle', 'none');
%         otherwise
%             agents(i).figure_handle_interests = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%     end
% end
% % draw target and map
% for i = 1:size(targets,1)
%     targets(i).figure_handle_interests = birds(targets(i).x, targets(i).y, targets(i).z, ...
%                 'body', color_pallete(i), 'yaw', rad2deg(targets(i).heading), 'scale', target_scale);
% end
% interests_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.inter_global, 'facealpha',1,'linestyle','none');
% oldcolormap = colormap('gray');
% colormap(flipud(oldcolormap));
% caxis([0 1]);
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;
% axis([-100 world.size_x+100 -100 world.size_y+100]);

% =======================================================================================================
% =======================================================================================================
% =======================================================================================================
%% debug figure
h2                  = figure('Name', 'Debug figure', 'NumberTitle', 'off');
h2.Position			= [1,41,1920,962]; % Biggest screen possible on my desktop size
drawBoundary;
hold on;

% CREATE legend dummy plots for legend entries on debug window
num_leg_entries = 5+numTargets;
hleglines = plot(1,ones(1,num_leg_entries));
UAV_state_string = strcat('\color[rgb]{',num2str(myBlue), '} CHASE \color{black}\\','\color[rgb]{',num2str(myDarkGreen), '} HOME \color{black}\\','\color[rgb]{',num2str(myLightGreen), '} RETURN');
leg_labels = {'Time: ',strcat("No.of UAVs: ",num2str(numAgentsArray(p))),strcat("Tether: ",num2str(maxHorizontalMovement),'m'),'UAV States:',UAV_state_string};
for rep = 1:numTargets
	% append extra legend info per target added
	leg_labels(end+1) = {strcat('\color[rgb]{',num2str(color_pallete(rep,:)), '} Bird flock'," ", num2str(rep), " energy: ")};
end

[hleg, icons] = legend(hleglines, leg_labels,'Location','northwest','FontSize',13);

set(findobj(icons,'Type','line','-not','Tag','dummy'),'Vis','off'); % finds any line objects which do not have the tag of 'dummy' and turns their visibility off
pos = get(findobj(icons,'Type','text'),'Pos'); % extracts position coordinates of legend entries
pos = cell2mat(pos);
h_txts = findobj(icons,'Type','text');

for indx = 1:length(pos(:,1))
% 	updated_pos{i} = [0.1 pos(i,2:3)];
	if indx == 5 
		continue
	end
	set(h_txts(indx),'Pos',[0.1 pos(indx,2:3)]);
end


% draw aircraft and map
% for i = 1:size(agents,1)
%     switch agents(i).type
%         case 1
%             agents(i).figure_handle_debug = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete_agent(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%         case 2
%             agents(i).figure_handle_debug = quadrotor(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', multirotor_scale, 'body', color_pallete_agent(i),...
%                      'boom', color_pallete(i), 'prop', color_pallete(i), 'linestyle', 'none');
%         otherwise
%             agents(i).figure_handle_debug = uav(agents(i).x, agents(i).y, agents(i).z,...
%                      'roll', rad2deg(agents(i).roll), 'pitch', rad2deg(agents(i).pitch),...
%                      'yaw', rad2deg(agents(i).yaw), 'scale', aircraft_scale, 'color', color_pallete_agent(i),...
%                      'wing', color_pallete(i), 'linestyle', 'none');
%     end
%     % draw position setpoint
%     agents(i).figure_handle_waypoint = plot3(agents(i).x, agents(i).y, agents(i).z, ...
%                 'o', 'markersize', 10, 'color', color_pallete_agent(i));
% end
% draw target and map
% for i = 1:size(targets,1)
%     targets(i).figure_handle_debug = birds(targets(i).x, targets(i).y, targets(i).z, ...
%                 'body', color_pallete(i), 'yaw', rad2deg(targets(i).heading), 'scale', target_scale);
%     draw position setpoint
%     targets(i).figure_handle_waypoint = plot3(targets(i).x, targets(i).y, targets(i).z, ...
%                 '^', 'markersize', 10, 'color', color_pallete(i));
% end
debug_map_handle = surf(world.meshgrid_XX,world.meshgrid_YY,world.prob_global(:,:,1), 'facealpha',1,'linestyle','none');
set(gcf,'color','w');
% Control the view aspect/angle of the debug figure:
view(viewAspectDesired);
% view([1,1,1]);
% an = annotation('textbox',[0.4 .85 0 0],'String',['Target energy reserve: ',num2str(targets(1).energy(1)),'Units'],'FontSize',12,'BackgroundColor','w','FitBoxToText','on');

states_y_pos = [0.85, 0.81, 0.77, 0.73];
state_options = ["Home", "Chase", "Standby", "Return"];
dodgy_yellow = [0.93, 0.86, 0.22];
an_colors = [myDarkGreen; mydarkblue; dodgy_yellow; myLightGreen];

% an_states = [];
% 
% for i = 1:length(states_y_pos)
% 	an_temp = annotation('textbox',[0.65 states_y_pos(i) 0 0],'String',state_options(i),'FontSize',12,'Color',an_colors(i,:),'BackgroundColor','w','FitBoxToText','on');
% 	an_states = [an_states; an_temp];
% end

oldcolormap = colormap('gray');
colormap(flipud(oldcolormap));
caxis([0 1]);
xlabel('x (m)','FontSize',20,'FontName','Arial','FontWeight','Bold');
ylabel('y (m)','FontSize',20,'FontName','Arial','FontWeight','Bold');
zlabel('z (m)','FontSize',20,'FontName','Arial','FontWeight','Bold');

daspect([2 2 1.7]); % relative sizing of scales in [x,y,z]
axis([0 world.size_x 0 world.size_y 0 20])





%%% COPY WORLD into subviews (for visual purposes)
% h3 = figure('Name', 'Debug figure', 'NumberTitle', 'off');
% hsub1 = subplot(2,1,1);
% hsub2 = subplot(2,1,2);
% copyobj(debug_map_handle, hsub1);
% view([-1,-1,0.55]);
% copyobj(debug_map_handle, hsub2);
% view([0, 0 , 1]);


