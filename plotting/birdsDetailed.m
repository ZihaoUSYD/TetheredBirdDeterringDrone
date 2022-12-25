function [h] = birdsDetailed(varargin)

bodyColor           = 0.5*[1 1 1];
boomColor           = 0.5*[1 1 1];
propColor           = 0.5*[1 1 1];
edgeColor           = 'k';
linestyle           = '-';
scale               = 1;
setroll             = false;
setpitch            = false;
setyaw              = false;

%% Parse Inputs
if nargin>2 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    assert(isscalar(x)==1, 'quad input x must be a scalar.')
    assert(isscalar(y)==1, 'quad input y must be a scalar.')
    assert(isscalar(z)==1, 'quad input z must be a scalar.')
else
    x = 0;
    y = 0;
    z = 0;
end

tmp = strncmpi(varargin, 'body', 4);
if any(tmp)
    bodyColor = varargin{find(tmp)+1};
end

tmp = strncmpi(varargin, 'lines', 5);
if any(tmp)
    linestyle = varargin{find(tmp)+1};
end

tmp = strncmpi(varargin,'linecol',7)|strncmpi(varargin,'edgecol',3); 
if any(tmp)
    edgeColor = varargin{find(tmp)+1}; 
end

tmp = strcmpi(varargin,'scale');  
if any(tmp) 
    scale = varargin{find(tmp)+1}; 
    assert(isscalar(scale)==1,'It may seem redundant, but the scale must be a scalar.')
end

tmp = strncmpi(varargin, 'yaw', 3);
if any(tmp)
    setyaw = true;
    yaw = varargin{find(tmp)+1};
end

tmp = strncmpi(varargin, 'moveAllow',9);
if any(tmp)
	move_flag = varargin{find(tmp)+1};
end

%% Dimensions
bodyWidth           = 0.1*scale;
bodyHeight          = 0.05*bodyWidth;
wingPart1Width		= 0.6*bodyWidth;
wingPart1Length		= 0.5*wingPart1Width;
wingPart1Disp		= 0.2*wingPart1Width;
wingPart2Width		= wingPart1Width;
wingPart2Length		= 1.1*wingPart2Width;
wingPart2Tip		= 1*wingPart1Width;
wingHeight			= 0.5*bodyHeight;
anhedralHeight		= 0.3*wingPart1Length;
% boomLength          = 0.25*scale;
% boomRadius          = 0.02*boomLength;
% propRadius          = 0.10*scale;
% propHeight          = 0.1*propRadius;

%% Determine if a figure is already open
initialHoldState    = 0;
SetView = isempty(get(gcf,'CurrentAxes'));

if ~SetView
    if ishold
        initialHoldState = 1;
    else
        cla;
        SetView     = true;
    end
end

% Toggle whether we want detailed bird models or not (simple models allow
% for the code to run faster/more efficiently)
detailsDesired = true;

%% Define scaling and relative positioning
% Body ratio scaling
w_s = 5; % width scaling of main body of bird (higher means thinner)
h_s = 3; % height scaling (higher means shorter)

% Add random noise/variation to bird position w.r.t to central bird
if move_flag
	x_upper = bodyWidth/3;
	x_lower = -bodyWidth/3;
	y_upper = bodyWidth/3;
	y_lower = -bodyWidth/3;
	z_upper	= 0.8; % upper bound of variation in metres
	z_lower	= -0.8; % lower bound of variation from central bird z in metres

else
	% no noise/fluttering if within range of given setpoint
	x_upper = 0;
	x_lower = 0;
	y_upper = 0;
	y_lower = 0;
	z_upper	= 0; % upper bound of variation in metres
	z_lower	= 0; % lower bound of variation from central bird z in metres
end
	x_noises = NaN(4,1);
	y_noises = NaN(4,1);
	z_noises = NaN(4,1);


for i = 1:length(z_noises)
	x_noises(i,1) = x_lower + (x_upper - x_lower)*rand(); % x displacement of bird i w.r.t central bird
	y_noises(i,1) = y_lower + (y_upper - y_lower)*rand(); % y displacement of bird i w.r.t central bird
	z_noises(i,1) = z_lower + (z_upper - z_lower)*rand(); % z displacement of bird i w.r.t central bird
end

%% Draw surfaces
% body
p1 = [x-bodyWidth/w_s, y-bodyWidth/h_s, z];
p2 = [x+bodyWidth/w_s, y-bodyWidth/h_s, z];
p3 = [x, y+bodyWidth, z];
p = [p1;p2;p3];
if ~initialHoldState
    hold on
end

h(21) = fill3(p(:,1),p(:,2),p(:,3)-bodyHeight,bodyColor);
h(22) = fill3(p(:,1),p(:,2),p(:,3)+bodyHeight,bodyColor);
h(23) = fill3([p1(1),p2(1),p2(1),p1(1)],...
             [p1(2),p2(2),p2(2),p1(2)],...
             [p1(3)+bodyHeight,p2(3)+bodyHeight,p2(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(24) = fill3([p1(1),p3(1),p3(1),p1(1)],...
             [p1(2),p3(2),p3(2),p1(2)],...
             [p1(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(25) = fill3([p2(1),p3(1),p3(1),p2(1)],...
             [p2(2),p3(2),p3(2),p2(2)],...
             [p2(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p2(3)-bodyHeight],...
             bodyColor);
		 
		 % BODY 1 - Add wing details of 
		 if detailsDesired
			 % wings of central bird (pw1_1_left = Coordinates of wing part 1 point 1
			 % on left side)
			 signs = [1, -1]; % correlate to left/right wings
			 for j = 1:length(signs)
				 pw11 = [x, y - 0.25*wingPart1Width, z];
				 pw12 = [x, y + 0.75*wingPart1Width, z];
				 pw13 = [x - signs(j)*wingPart1Length, y - 0.25*wingPart1Width + wingPart1Disp, z + anhedralHeight];
				 pw14 = [x - signs(j)*wingPart1Length, y + 0.75*wingPart1Width + wingPart1Disp, z + anhedralHeight];
				 pw2_tip  = [x - signs(j)*wingPart1Length - signs(j)*wingPart2Length, y + wingPart1Disp - wingPart2Tip, z - anhedralHeight];
				 
				 pw1_coords = [pw11; pw12; pw14; pw13]; % wing part 1 rhombus (left side)
				 pw2_coords = [pw13; pw14; pw2_tip]; % wing part 2 trianlge (left side)
				 rep = 11; % This number just represents the number of shapes being mirrored across the left and right sides --> just for plotting handle purposes
				 
				 % Rhomboids top & bottom
				 h(26+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)-wingHeight, bodyColor); % bottom rhombus - left side
				 h(27+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)+wingHeight, bodyColor); % top rhombus - left side
				 
				 % Vertical walls of rhomboids
				 h(28+rep*(j-1)) = fill3([pw11(1),pw12(1),pw12(1),pw11(1)],...
					 [pw11(2),pw12(2),pw12(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw12(3)+wingHeight,pw12(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 h(29+rep*(j-1)) = fill3([pw12(1),pw14(1),pw14(1),pw12(1)],...
					 [pw12(2),pw14(2),pw14(2),pw12(2)],...
					 [pw12(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw12(3)-wingHeight],...
					 bodyColor);
				 h(30+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(31+rep*(j-1)) = fill3([pw11(1),pw13(1),pw13(1),pw11(1)],...
					 [pw11(2),pw13(2),pw13(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 % Wing tips
				 h(32+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)-wingHeight, bodyColor); % bottom rhombus
				 h(33+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)+wingHeight, bodyColor); % top rhombus
				 h(34+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(35+rep*(j-1)) = fill3([pw14(1),pw2_tip(1),pw2_tip(1),pw14(1)],...
					 [pw14(2),pw2_tip(2),pw2_tip(2),pw14(2)],...
					 [pw14(3)+wingHeight,pw2_tip(3)+wingHeight,pw2_tip(3)-wingHeight,pw14(3)-wingHeight],...
					 bodyColor);
				 h(36+rep*(j-1)) = fill3([pw2_tip(1),pw13(1),pw13(1),pw2_tip(1)],...
					 [pw2_tip(2),pw13(2),pw13(2),pw2_tip(2)],...
					 [pw2_tip(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw2_tip(3)-wingHeight],...
					 bodyColor);
			 end
			 
		 end

		 
% body 2
bodyNoise2 = [x_noises(1), 2*bodyWidth + y_noises(1), z_noises(1)];
p1 = [x-bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise2;
p2 = [x+bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise2;
p3 = [x, y+bodyWidth, z] + bodyNoise2;
p = [p1;p2;p3];
if ~initialHoldState
    hold on
end

h(1) = fill3(p(:,1),p(:,2),p(:,3)-bodyHeight,bodyColor);
h(2) = fill3(p(:,1),p(:,2),p(:,3)+bodyHeight,bodyColor);
h(3) = fill3([p1(1),p2(1),p2(1),p1(1)],...
             [p1(2),p2(2),p2(2),p1(2)],...
             [p1(3)+bodyHeight,p2(3)+bodyHeight,p2(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(4) = fill3([p1(1),p3(1),p3(1),p1(1)],...
             [p1(2),p3(2),p3(2),p1(2)],...
             [p1(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(5) = fill3([p2(1),p3(1),p3(1),p2(1)],...
             [p2(2),p3(2),p3(2),p2(2)],...
             [p2(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p2(3)-bodyHeight],...
             bodyColor);
		 
		 % BODY 2 - Add wing details (start from fig handle #48)
		 if detailsDesired
			 % wings of central bird (pw1_1_left = Coordinates of wing part 1 point 1
			 % on left side)
			 signs = [1, -1]; % correlate to left/right wings
			 for j = 1:length(signs)
				 pw11 = [x, y - 0.25*wingPart1Width, z]+bodyNoise2;
				 pw12 = [x, y + 0.75*wingPart1Width, z]+bodyNoise2;
				 pw13 = [x - signs(j)*wingPart1Length, y - 0.25*wingPart1Width + wingPart1Disp, z + anhedralHeight]+bodyNoise2;
				 pw14 = [x - signs(j)*wingPart1Length, y + 0.75*wingPart1Width + wingPart1Disp, z + anhedralHeight]+bodyNoise2;
				 pw2_tip  = [x - signs(j)*wingPart1Length - signs(j)*wingPart2Length, y + wingPart1Disp - wingPart2Tip, z - anhedralHeight]+ bodyNoise2;
				 
				 pw1_coords = [pw11; pw12; pw14; pw13]; % wing part 1 rhombus (left side)
				 pw2_coords = [pw13; pw14; pw2_tip]; % wing part 2 trianlge (left side)
				 rep = 11; % This number just represents the number of shapes being mirrored across the left and right sides --> just for plotting handle purposes
				 
				 % Rhomboids top & bottom
				 h(48+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)-wingHeight, bodyColor); % bottom rhombus - left side
				 h(49+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)+wingHeight, bodyColor); % top rhombus - left side
				 
				 % Vertical walls of rhomboids
				 h(50+rep*(j-1)) = fill3([pw11(1),pw12(1),pw12(1),pw11(1)],...
					 [pw11(2),pw12(2),pw12(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw12(3)+wingHeight,pw12(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 h(51+rep*(j-1)) = fill3([pw12(1),pw14(1),pw14(1),pw12(1)],...
					 [pw12(2),pw14(2),pw14(2),pw12(2)],...
					 [pw12(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw12(3)-wingHeight],...
					 bodyColor);
				 h(52+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(53+rep*(j-1)) = fill3([pw11(1),pw13(1),pw13(1),pw11(1)],...
					 [pw11(2),pw13(2),pw13(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 % Wing tips
				 h(54+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)-wingHeight, bodyColor); % bottom rhombus
				 h(55+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)+wingHeight, bodyColor); % top rhombus
				 h(56+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(57+rep*(j-1)) = fill3([pw14(1),pw2_tip(1),pw2_tip(1),pw14(1)],...
					 [pw14(2),pw2_tip(2),pw2_tip(2),pw14(2)],...
					 [pw14(3)+wingHeight,pw2_tip(3)+wingHeight,pw2_tip(3)-wingHeight,pw14(3)-wingHeight],...
					 bodyColor);
				 h(58+rep*(j-1)) = fill3([pw2_tip(1),pw13(1),pw13(1),pw2_tip(1)],...
					 [pw2_tip(2),pw13(2),pw13(2),pw2_tip(2)],...
					 [pw2_tip(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw2_tip(3)-wingHeight],...
					 bodyColor);
			 end
			 
		 end
		 
		 
         
% body 3
bodyNoise3 = [x_noises(2), -2*bodyWidth + y_noises(2),z_noises(2)];
p1 = [x-bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise3;
p2 = [x+bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise3;
p3 = [x, y+bodyWidth, z] + bodyNoise3;
p = [p1;p2;p3];
if ~initialHoldState
    hold on
end

h(6) = fill3(p(:,1),p(:,2),p(:,3)-bodyHeight,bodyColor);
h(7) = fill3(p(:,1),p(:,2),p(:,3)+bodyHeight,bodyColor);
h(8) = fill3([p1(1),p2(1),p2(1),p1(1)],...
             [p1(2),p2(2),p2(2),p1(2)],...
             [p1(3)+bodyHeight,p2(3)+bodyHeight,p2(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(9) = fill3([p1(1),p3(1),p3(1),p1(1)],...
             [p1(2),p3(2),p3(2),p1(2)],...
             [p1(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(10) = fill3([p2(1),p3(1),p3(1),p2(1)],...
             [p2(2),p3(2),p3(2),p2(2)],...
             [p2(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p2(3)-bodyHeight],...
             bodyColor);
		 
		 % BODY 3 - Add wing details (start from fig handle #69)
		 if detailsDesired
			 % wings of central bird (pw1_1_left = Coordinates of wing part 1 point 1
			 % on left side)
			 signs = [1, -1]; % correlate to left/right wings
			 for j = 1:length(signs)
				 pw11 = [x, y - 0.25*wingPart1Width, z] + bodyNoise3;
				 pw12 = [x, y + 0.75*wingPart1Width, z] + bodyNoise3;
				 pw13 = [x - signs(j)*wingPart1Length, y - 0.25*wingPart1Width + wingPart1Disp, z + anhedralHeight] + bodyNoise3;
				 pw14 = [x - signs(j)*wingPart1Length, y + 0.75*wingPart1Width + wingPart1Disp, z + anhedralHeight] + bodyNoise3;
				 pw2_tip  = [x - signs(j)*wingPart1Length - signs(j)*wingPart2Length, y + wingPart1Disp - wingPart2Tip, z - anhedralHeight] + bodyNoise3;
				 
				 pw1_coords = [pw11; pw12; pw14; pw13]; % wing part 1 rhombus (left side)
				 pw2_coords = [pw13; pw14; pw2_tip]; % wing part 2 trianlge (left side)
				 rep = 11; % This number just represents the number of shapes being mirrored across the left and right sides --> just for plotting handle purposes
				 
				 % Rhomboids top & bottom
				 h(70+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)-wingHeight, bodyColor); % bottom rhombus - left side
				 h(71+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)+wingHeight, bodyColor); % top rhombus - left side
				 
				 % Vertical walls of rhomboids
				 h(72+rep*(j-1)) = fill3([pw11(1),pw12(1),pw12(1),pw11(1)],...
					 [pw11(2),pw12(2),pw12(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw12(3)+wingHeight,pw12(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 h(73+rep*(j-1)) = fill3([pw12(1),pw14(1),pw14(1),pw12(1)],...
					 [pw12(2),pw14(2),pw14(2),pw12(2)],...
					 [pw12(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw12(3)-wingHeight],...
					 bodyColor);
				 h(74+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(75+rep*(j-1)) = fill3([pw11(1),pw13(1),pw13(1),pw11(1)],...
					 [pw11(2),pw13(2),pw13(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 % Wing tips
				 h(76+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)-wingHeight, bodyColor); % bottom rhombus
				 h(77+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)+wingHeight, bodyColor); % top rhombus
				 h(78+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(79+rep*(j-1)) = fill3([pw14(1),pw2_tip(1),pw2_tip(1),pw14(1)],...
					 [pw14(2),pw2_tip(2),pw2_tip(2),pw14(2)],...
					 [pw14(3)+wingHeight,pw2_tip(3)+wingHeight,pw2_tip(3)-wingHeight,pw14(3)-wingHeight],...
					 bodyColor);
				 h(80+rep*(j-1)) = fill3([pw2_tip(1),pw13(1),pw13(1),pw2_tip(1)],...
					 [pw2_tip(2),pw13(2),pw13(2),pw2_tip(2)],...
					 [pw2_tip(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw2_tip(3)-wingHeight],...
					 bodyColor);
			 end
			 
		 end
		 
		 
         
% body 4
bodyNoise4 = [2*bodyWidth + x_noises(3), y_noises(3),z_noises(3)];
p1 = [x-bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise4;
p2 = [x+bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise4;
p3 = [x, y+bodyWidth, z] + bodyNoise4;
p = [p1;p2;p3];
if ~initialHoldState
    hold on
end

h(11) = fill3(p(:,1),p(:,2),p(:,3)-bodyHeight,bodyColor);
h(12) = fill3(p(:,1),p(:,2),p(:,3)+bodyHeight,bodyColor);
h(13) = fill3([p1(1),p2(1),p2(1),p1(1)],...
             [p1(2),p2(2),p2(2),p1(2)],...
             [p1(3)+bodyHeight,p2(3)+bodyHeight,p2(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(14) = fill3([p1(1),p3(1),p3(1),p1(1)],...
             [p1(2),p3(2),p3(2),p1(2)],...
             [p1(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(15) = fill3([p2(1),p3(1),p3(1),p2(1)],...
             [p2(2),p3(2),p3(2),p2(2)],...
             [p2(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p2(3)-bodyHeight],...
             bodyColor);
		 
		 % BODY 4 - Add wing details (start from fig handle #69)
		 if detailsDesired
			 % wings of central bird (pw1_1_left = Coordinates of wing part 1 point 1
			 % on left side)
			 signs = [1, -1]; % correlate to left/right wings
			 for j = 1:length(signs)
				 pw11 = [x, y - 0.25*wingPart1Width, z] + bodyNoise4;
				 pw12 = [x, y + 0.75*wingPart1Width, z] + bodyNoise4;
				 pw13 = [x - signs(j)*wingPart1Length, y - 0.25*wingPart1Width + wingPart1Disp, z + anhedralHeight] + bodyNoise4;
				 pw14 = [x - signs(j)*wingPart1Length, y + 0.75*wingPart1Width + wingPart1Disp, z + anhedralHeight] + + bodyNoise4;
				 pw2_tip  = [x - signs(j)*wingPart1Length - signs(j)*wingPart2Length, y + wingPart1Disp - wingPart2Tip, z - anhedralHeight] + bodyNoise4;
				 
				 pw1_coords = [pw11; pw12; pw14; pw13]; % wing part 1 rhombus (left side)
				 pw2_coords = [pw13; pw14; pw2_tip]; % wing part 2 trianlge (left side)
				 rep = 11; % This number just represents the number of shapes being mirrored across the left and right sides --> just for plotting handle purposes
				 
				 % Rhomboids top & bottom
				 h(92+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)-wingHeight, bodyColor); % bottom rhombus - left side
				 h(93+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)+wingHeight, bodyColor); % top rhombus - left side
				 
				 % Vertical walls of rhomboids
				 h(94+rep*(j-1)) = fill3([pw11(1),pw12(1),pw12(1),pw11(1)],...
					 [pw11(2),pw12(2),pw12(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw12(3)+wingHeight,pw12(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 h(95+rep*(j-1)) = fill3([pw12(1),pw14(1),pw14(1),pw12(1)],...
					 [pw12(2),pw14(2),pw14(2),pw12(2)],...
					 [pw12(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw12(3)-wingHeight],...
					 bodyColor);
				 h(96+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(97+rep*(j-1)) = fill3([pw11(1),pw13(1),pw13(1),pw11(1)],...
					 [pw11(2),pw13(2),pw13(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 % Wing tips
				 h(98+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)-wingHeight, bodyColor); % bottom rhombus
				 h(99+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)+wingHeight, bodyColor); % top rhombus
				 h(100+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(101+rep*(j-1)) = fill3([pw14(1),pw2_tip(1),pw2_tip(1),pw14(1)],...
					 [pw14(2),pw2_tip(2),pw2_tip(2),pw14(2)],...
					 [pw14(3)+wingHeight,pw2_tip(3)+wingHeight,pw2_tip(3)-wingHeight,pw14(3)-wingHeight],...
					 bodyColor);
				 h(102+rep*(j-1)) = fill3([pw2_tip(1),pw13(1),pw13(1),pw2_tip(1)],...
					 [pw2_tip(2),pw13(2),pw13(2),pw2_tip(2)],...
					 [pw2_tip(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw2_tip(3)-wingHeight],...
					 bodyColor);
			 end
			 
		 end
		 
		 	 
         
% body 5
bodyNoise5 = [-2*bodyWidth + x_noises(4), y_noises(4),  z_noises(4)];
p1 = [x-bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise5;
p2 = [x+bodyWidth/w_s, y-bodyWidth/h_s, z] + bodyNoise5;
p3 = [x, y+bodyWidth, z] + bodyNoise5;
p = [p1;p2;p3];
if ~initialHoldState
    hold on
end

h(16) = fill3(p(:,1),p(:,2),p(:,3)-bodyHeight,bodyColor);
h(17) = fill3(p(:,1),p(:,2),p(:,3)+bodyHeight,bodyColor);
h(18) = fill3([p1(1),p2(1),p2(1),p1(1)],...
             [p1(2),p2(2),p2(2),p1(2)],...
             [p1(3)+bodyHeight,p2(3)+bodyHeight,p2(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(19) = fill3([p1(1),p3(1),p3(1),p1(1)],...
             [p1(2),p3(2),p3(2),p1(2)],...
             [p1(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p1(3)-bodyHeight],...
             bodyColor);
h(20) = fill3([p2(1),p3(1),p3(1),p2(1)],...
             [p2(2),p3(2),p3(2),p2(2)],...
             [p2(3)+bodyHeight,p3(3)+bodyHeight,p3(3)-bodyHeight,p2(3)-bodyHeight],...
             bodyColor);
		 
		 % BODY 5 - Add wing details (start from fig handle #69)
		 if detailsDesired
			 % wings of central bird (pw1_1_left = Coordinates of wing part 1 point 1
			 % on left side)
			 signs = [1, -1]; % correlate to left/right wings
			 for j = 1:length(signs)
				 pw11 = [x, y - 0.25*wingPart1Width, z] + bodyNoise5;
				 pw12 = [x, y + 0.75*wingPart1Width, z] + + bodyNoise5;
				 pw13 = [x - signs(j)*wingPart1Length, y - 0.25*wingPart1Width + wingPart1Disp, z + anhedralHeight] + bodyNoise5;
				 pw14 = [x - signs(j)*wingPart1Length, y + 0.75*wingPart1Width + wingPart1Disp, z + anhedralHeight] + + bodyNoise5;
				 pw2_tip  = [x - signs(j)*wingPart1Length - signs(j)*wingPart2Length, y + wingPart1Disp - wingPart2Tip, z - anhedralHeight] + bodyNoise5;
				 
				 pw1_coords = [pw11; pw12; pw14; pw13]; % wing part 1 rhombus (left side)
				 pw2_coords = [pw13; pw14; pw2_tip]; % wing part 2 trianlge (left side)
				 rep = 11; % This number just represents the number of shapes being mirrored across the left and right sides --> just for plotting handle purposes
				 
				 % Rhomboids top & bottom
				 h(114+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)-wingHeight, bodyColor); % bottom rhombus - left side
				 h(115+rep*(j-1)) = fill3(pw1_coords(:,1),pw1_coords(:,2),pw1_coords(:,3)+wingHeight, bodyColor); % top rhombus - left side
				 
				 % Vertical walls of rhomboids
				 h(116+rep*(j-1)) = fill3([pw11(1),pw12(1),pw12(1),pw11(1)],...
					 [pw11(2),pw12(2),pw12(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw12(3)+wingHeight,pw12(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 h(117+rep*(j-1)) = fill3([pw12(1),pw14(1),pw14(1),pw12(1)],...
					 [pw12(2),pw14(2),pw14(2),pw12(2)],...
					 [pw12(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw12(3)-wingHeight],...
					 bodyColor);
				 h(118+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(119+rep*(j-1)) = fill3([pw11(1),pw13(1),pw13(1),pw11(1)],...
					 [pw11(2),pw13(2),pw13(2),pw11(2)],...
					 [pw11(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw11(3)-wingHeight],...
					 bodyColor);
				 % Wing tips
				 h(120+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)-wingHeight, bodyColor); % bottom rhombus
				 h(121+rep*(j-1)) = fill3(pw2_coords(:,1),pw2_coords(:,2),pw2_coords(:,3)+wingHeight, bodyColor); % top rhombus
				 h(122+rep*(j-1)) = fill3([pw13(1),pw14(1),pw14(1),pw13(1)],...
					 [pw13(2),pw14(2),pw14(2),pw13(2)],...
					 [pw13(3)+wingHeight,pw14(3)+wingHeight,pw14(3)-wingHeight,pw13(3)-wingHeight],...
					 bodyColor);
				 h(123+rep*(j-1)) = fill3([pw14(1),pw2_tip(1),pw2_tip(1),pw14(1)],...
					 [pw14(2),pw2_tip(2),pw2_tip(2),pw14(2)],...
					 [pw14(3)+wingHeight,pw2_tip(3)+wingHeight,pw2_tip(3)-wingHeight,pw14(3)-wingHeight],...
					 bodyColor);
				 h(124+rep*(j-1)) = fill3([pw2_tip(1),pw13(1),pw13(1),pw2_tip(1)],...
					 [pw2_tip(2),pw13(2),pw13(2),pw2_tip(2)],...
					 [pw2_tip(3)+wingHeight,pw13(3)+wingHeight,pw13(3)-wingHeight,pw2_tip(3)-wingHeight],...
					 bodyColor);
			 end
			 
		 end
		 
		 
		 
		 
%% Set roll, pitch, yaw
if setyaw
    rotate(h,[0 0 1],yaw-90,[x y z]);
end

%% Set view
% if SetView
%     view([140 30]);
%     axis tight equal
%     lighting gouraud
%     camlight
% else
%     axis auto
% end
