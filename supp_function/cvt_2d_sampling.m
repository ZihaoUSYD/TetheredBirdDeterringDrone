function [g] = cvt_2d_sampling ( g_num, it_num, s_num, height, width, sim_num)

%*****************************************************************************80
%
%% CVT_2D_SAMPLING carries out the Lloyd algorithm in a 2D unit box.
%
%  Discussion:
%
%    This program is a variation of the CVT_2D_LLOYD method.
%
%    Instead of using an exact technique to determine the Voronoi
%    regions, it uses sampling.
%
%    MATLAB highly inconvenienced me (let me put it nicely) by removing
%    the DSEARCH function, making all my CVT codes fail until I found
%    the replacement function.  It would have killed them to leave the
%    old thing available, so my codes would still run?
%    JVB, 06 June 2016.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    13/5/2021 by Joshua Trethowan
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    integer G_NUM, the number of generators.
%    A value of 50 is reasonable.
%
%    integer IT_NUM, the number of CVT iterations.
%    A value of 20 or 50 might be reasonable.
%
%    integer S_NUM, the number of sample points to use
%    when estimating the Voronoi regions.
%    A value of 1,000 is too low.  A value of 1,000,000 is somewhat high.
%
%   fprintf ( 1, '\n' );
%   fprintf ( 1, 'CVT_2D_SAMPLING\n' );
%   fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
  fprintf ( 1, '  Use sampling to approximate Lloyd''s algorithm\n' );
  fprintf ( 1, '  in the 2D %d m x %d m rectangular grid.\n',width, height);

  if ( nargin < 1 )
    g_num = input ( '  Enter number of generators: ' );
  elseif ( ischar ( g_num ) )
    g_num = str2num ( g_num );
  end

  if ( nargin < 2 ) 
    it_num = input ( '  Enter number of iterations: ' );
  elseif ( ischar ( it_num ) )
    it_num = str2num ( it_num );
  end

  if ( nargin < 3 ) 
    s_num = input ( '  Enter number of sample points: ' );
  elseif ( ischar ( s_num ) )
    s_num = str2num ( s_num );
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, '  Number of generators is %d\n', g_num );
  fprintf ( 1, '  Number of iterations is %d\n', it_num );
  fprintf ( 1, '  Number of samples is %d\n', s_num );
%

% Introduce scaling factors
% height = 400;
% width = 700;
%  Initialize the generators.
bound_lower = 0; % fraction of height or width, was 0.45
bound_upper = 1; % fraction of height or width, was 0.55
  g_height = bound_lower*height + (bound_upper*height - bound_lower*height)*rand ( g_num, 1 );
  g_width = bound_lower*width + (bound_upper*width - bound_lower*width)*rand (g_num, 1);
  
%   g_height = linspace(0, height, g_num)';
%   g_width = linspace(0, width, g_num)';
  
  % Add some noise to the data
%   g_height = g_height + rand(length(g_height),1)*0.1*height;
%   g_width = g_width + rand(length(g_width),1)*0.1*width;
  

g = [g_width, g_height];
%
%  Carry out the iteration.
%
  step = 1 : it_num;
  e = nan ( it_num, 1 );
  gm = nan ( it_num, 1 );

  f = figure ('Visible','off'); % set visibility to 'on' if you want to see how the points conform dynamically
  
  for it = 1 : it_num
%
%  Compute the Delaunay triangle information T for the current nodes.
%
pointsAreCollinear = @(xy) rank(xy(2:end,:) - xy(1,:)) == 1;
	if pointsAreCollinear(g) == true
		% already collinear, can return early
		return
	else
		t = delaunay ( g(:,1), g(:,2) );
	end
%
%  Display the Voronoi cells.
%
%     voronoi ( g(:,1), g(:,2), t );
% %     title_string = sprintf ( 'Voronoi Diagram: Agent Grid Coordinates, step %d', it );
% %     title ( title_string );
%     axis equal
%     axis ( [  0.0, width, 0.0, height ] )
%     axis square
% 	xlabel('x (m)');
% 	ylabel('y (m)');
%     drawnow
%  GraphGood();
% 	set(f, 'visible', 'on','Position',[100 100 1600 1600]);
%
%  Display the Delaunay triangulation.
%
%     figure ( 2 );
%     trimesh ( t, g(:,1), g(:,2), zeros(g_num,1), 'EdgeColor', 'r' )
%     hold on
%     scatter ( g(:,1), g(:,2), 'k.' )
%     title_string = sprintf ( 'Delaunay, step %d', it );
%     title ( title_string );
% %     axis ( [  0.0, width, 0.0, height ] )
%     axis square
%     view ( 2 )
%     hold off
%
%  Generate sample points.
%  New option for fixed grid sampling.
%
	if ( false )
		s = rand ( s_num, 2 );
	else
% 		s2 = max ( s2, 2 );
		%       s_num = s2 * s2;
		if height == width
			% Property is a square
			s_width = floor ( sqrt ( s_num ) );
			s_height = s_width;
		else
			hw_ratio = height/width; % relative height to weight ratio will be used to scale where the samples should be taken
			s_width = floor(sqrt(1/hw_ratio*s_num));
			s_height = floor(s_width * hw_ratio);
		end
		
		s_num = s_height * s_width;
		x_divs = linspace(0, width, s_width);
		y_divs = linspace(0, height, s_height);
		[ sx, sy ] = meshgrid (x_divs, y_divs );
		sx = reshape ( sx, s_width * s_height, 1 );
		sy = reshape ( sy, s_width * s_height, 1 );
		s = [ sx, sy ];
	end
	%
%  For each sample point, find K, the index of the nearest generator.
%  We do this efficiently by using the Delaunay information with
%  Matlab's DSEARCH command, rather than a brute force nearest neighbor
%  computation.  Also, DSEARCH has been removed, so we need DSEARCHN.
%  
    k = dsearchn ( g, t, s );

    m(:,1) = accumarray ( k, ones(s_num,1) );

    g_new(:,1) = accumarray ( k, s(:,1) ) ./ m(:,1); % effectively keeps reducing the relative distance between each generator so they are as equally spaced as possible
    g_new(:,2) = accumarray ( k, s(:,2) ) ./ m(:,1);
%
%  Compute the energy.
%
    e(it,1) = sum ( ( s(:,1) - g(k(:,1),1) ).^2 ...
                  + ( s(:,2) - g(k(:,1),2) ).^2 ) / s_num;
%
%  Display the energy.
%
%     figure ( 3 )
%     plot ( step, log ( e ), 'm-*' )
%     title ( 'Log (Energy)' )
%     xlabel ( 'Step' )
%     ylabel ( 'Energy' )
%     grid
% %
%  Compute the generator motion.
%
    gm(it,1) = sum ( ( g_new(:,1) - g(:,1) ).^2 ...
                   + ( g_new(:,2) - g(:,2) ).^2 ) / g_num;
%
%  Display the generator motion.
%
%     figure ( 4 )
%     plot ( step, log ( gm ), 'm-*' )
%     title ( 'Log (Average generator motion)' )
%     xlabel ( 'Step' )
%     ylabel ( 'Energy' )
%     grid
%
%  Update the generators.
%
    g = g_new;

  end
%
%  Save final graphics.
%
%   figure(f) % uncomment this if you want to see Voronoi diagram pop up

% f = figure ('Visible','off'); % set visibility to 'on' if you want to see how the points conform dynamically
%     voronoi ( g(:,1), g(:,2), t );
%     title_string = sprintf ( 'Voronoi Diagram: Agent Grid Coordinates, step %d', it );
%     title ( title_string );
%     axis equal
%     axis ( [  0.0, width, 0.0, height ] )
%     axis square
% 
% filename = sprintf('voronoi_%d_agents_sim%d.png',g_num, sim_num);
% print ( '-dpng', filename );
% fprintf ( 1, '  Graphics saved as "%s"\n', filename );
% close(f);

%   figure ( 2 )
%   filename = 'delaunay.png';
%   print ( '-dpng', filename );
%   fprintf ( 1, '  Graphics saved as "%s"\n', filename );

%   figure ( 3 )
%   filename = 'energy.png';
%   print ( '-dpng', filename );
%   fprintf ( 1, '  Graphics saved as "%s"\n', filename );
% 
%   figure ( 4 )
%   filename = 'motion.png';
%   print ( '-dpng', filename );
%   fprintf ( 1, '  Graphics saved as "%s"\n', filename );
% %
%  Terminate.
%
%   fprintf ( 1, '\n' );
%   fprintf ( 1, 'CVT_2D_SAMPLING\n' );
%   fprintf ( 1, '  Normal end of execution.\n' );

  return
end

