function [] = GraphGood()
    %% Plot parameters
%     myred           = [216 30 49]/255;
%     myblue          = [0 72 255]/255;
%     myred           = [255 0 0]/255;
% 	mypink			= [255,20,147]/255;
%     myblue          = [0 0 255]/255;
%     myblack         = [0 0 0]/255;
% %     mygreen         = [0 200 0]/255;
%     mygreen         = [0 255 0]/255;
% 	mydarkgreen		= [46,139,87]/255;
%     mycyan          = [2 169 226]/255;
%     myyellow        = [249 190 0]/255;
%     mygray          = [89 89 89]/255;
% 	myorange		= [255,140,0]/255;
% 	myviolet		= [148,0,211]/255;
% 	mydarkred		= [139,0,0]/255;
% 	mytan			= [222,184,135]/255;
% 	myaqua			= [127,255,212]/255;
	
    alw             = 1;                        % AxesLineWidth
    fsz             = 12;                       % Fontsize
    lw              = 2;                        % LineWidth
    msz             = 15;                       % MarkerSize
    set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
    set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
    set(0,'defaultAxesLineWidth',alw);           % set the default line width to lw
    set(0,'defaultAxesFontSize',fsz);         % set the default line marker size to msz
    set(0,'defaulttextInterpreter','latex')
    set(0,'defaultAxesTickLabelInterpreter', 'latex');
    set(0,'defaultLegendInterpreter', 'latex');
    set(0,'defaultFigureColor','w');
    set(0,'defaultAxesColor','w');

    %% How to plot
    % hFig            = figure();
    % hold on; grid on; grid minor; axis tight; box on; 
    % 
    % plot(x,y)
    % 
    % xlabel('xlabel ($s$)')
    % ylabel('ylabel ($\phi$)')
    % 
    % set(gca,'GridLineStyle','-')
    % set(gca,'MinorGridLineStyle','-')
    % set(gca,'GridColor','k')
    % set(gca,'MinorGridColor','k')
    % 
    % hleg = legend('$\phi$','$\theta$','$\psi$');
    % set(hleg,'EdgeColor',hleg.Color);
    % set(hleg,'Location','best');
    % set(hleg,'Interpreter','latex')
end