function plotFigure(xData, yData, ColorPalette, xLabel, yLabel, titleString, legendString)

	plot(xData,yData,'Color',ColorPalette,'LineStyle','-','LineWidth',2);
	hold on
	xlabel(xLabel);
	ylabel(yLabel);
	legendString5 = sprintf("%d agent/s",numAgentsArray(p));
	legendStringArray5 = [legendStringArray5; legendString5];
	titleString5 = sprintf("%d m x %d m grid: Simulation %d", property_size(1),property_size(2),sim);
		
end