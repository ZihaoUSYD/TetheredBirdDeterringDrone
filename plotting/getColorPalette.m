function color_pallete_agent = getColorPalette(agent_state,flag1,flag2,flag3,color_pallete_agent_state1,color_pallete_agent_state2,color_pallete_agent_state3, color_pallete_agent_state4)
%% Get the correct colour palette depending on the agent state
if agent_state == flag1
	color_pallete_agent = color_pallete_agent_state1;
elseif agent_state == flag2
	color_pallete_agent = color_pallete_agent_state2;
elseif agent_state == flag3
	color_pallete_agent = color_pallete_agent_state4;
else
	color_pallete_agent = color_pallete_agent_state3;
end