classdef World
    % ===============================================================
    % Class - World
    % This class contains all the properties the simulation world
    %
    % Methods: Protected
    %           [obj] = getCentre(obj)
    %           [obj] = getNumberOfTiles(obj)
    %           [obj] = getMeshGrids(obj)
    %           [obj] = getGridCentres(obj)
    %           [obj] = initProbability(obj)
    %           [obj] = initUncertainty(obj)
    %           [obj] = initInterests(obj)
    %           [obj] = initRHC(obj)
    %          Private
    %           [obj] = expandArrays(obj, iterations)
    %           [obj] = updateGlobalInterests(obj, targets, cumulation)
    %           [obj] = updateGlobalUncertainty(obj, agents, cumulation)
    %           [obj] = updateGlobalProbability(obj, agents, cumulation)
    %
    % Created by: Zihao Wang
    % Last updated: 19-07-2018
    % ===============================================================
    properties
        % geometry (input)
        shape       % (rectangle|square|circle)
        grid_size   % size of the indivdual square grids
		buffer_size % [x,y] buffer size outside of property
		
        % rectangle
        size_x      % size of the world along x dimension
        size_y      % size of the world along y dimension
        % circle and square
        size_other  % edge length or diameter
        
        % geometry (calculated)
        world_centre
        number_of_tiles_x
        number_of_tiles_y
        number_of_tiles
        grid_centre_x
        grid_centre_y
		property_x_index % used for creating interest bound box
		property_y_index % used for creating interest bound box
        
        % meshgrids
        meshgrid_XX
        meshgrid_YY
        
        % probability (range 0-1)
        % 1: target is 100% there
        % 2: target is 100% not there
        prob0                % initial probability
        prob_nom             % nominal probabiltiy value (probabiltiy at infinity)
        prob_tau             % time constant - rate at which probability approach nominal value
        prob_detect          % maximum probability of grids in sensor and target detected
        prob_detect_tau      % time constant - rate at which probability approach maximum probability detection rate
        prob_nodetect        % probability of grids in sensor but no target detected
        prob_nodetect_tau    % time constant - rate at which probability approach minimum probability detection rate
        prob_range           % specify a circle around a potential target where probability will increase
        
        % uncertainty (range 0-1)
        % 1: no information about the world
        % 0: just observed the world
        uncert0     % initial uncertainty
        uncert_nom  % nominal uncertainty (uncertainty at infinity)
        uncert_tau  % time constant - rate at which uncertainty approach nominal value
        
        % interests (range 0-1)
        % 1: target will want to stay there
        % 0: target will want to avoid at all cost
        inter0              % initial interests
        inter_wild          % initial interests outside the vineyard
        inter_nom           % nominal interests
        inter_tau           % time constant - rate at which interests approach nominal value
        inter_nogo          % the interests in the nogo radius
        inter_nogo_radius   % the radius around the agent where the interests is instantly 0
        inter_alert_radius  % the radius around the agent where the interests will start to decrease
        
        % global knowledge
        prob_global     % global knowledge of probability
        uncert_global   % global knowledge of uncertainty
        inter_global    % global knowledge of interests
        
        % Receding Horizon Control
        RHC_steps
    end
    
    methods
        %% constructor
        function obj = World(shape,size,buffer_size,grid_size)
            if nargin > 0
                obj.shape       = shape;
                obj.grid_size   = grid_size;
				obj.buffer_size	= buffer_size;
				
                % case 1 rectangular property
                if strcmp(shape, 'rectangle')
                    obj.size_x  = size(1);
                    obj.size_y  = size(2);
                % case 2 square and circular property
                elseif strcmp(shape, 'square') || strcmp(shape, 'circle')
                    obj.size_other = size;
                    obj.size_x  = size;
                    obj.size_y  = size;
                else
                    warning('Invalid world shape argument, defaults to square');
                    obj.shape   = 'square';
                    obj.size_other = size(1);
                    obj.size_x  = size(1);
                    obj.size_y  = size(1);
                end
                
            else % default world shape is square with size 100, grid size 1
                obj.shape       = 'square';
				obj.buffer_size = [50, 50];
                obj.size_other  = 100;
                obj.size_x      = 100;
                obj.size_y      = 100;
                obj.grid_size   = 1;
            end
            
            % get calculated geometry
            obj = getCentre(obj);
            obj = getNumberOfTiles(obj);
            obj = getMeshGrids(obj);
            obj = getGridCentres(obj);
            
            % initialise target probability, environment uncertainty and interests maps
            % probability gives indication about where the potential target is
            % uncertainty indicates wether the information about an area is complete
            % uncertainty is 0 regardless of whether a target is present
            % as soon as it is within the UAV's FOV, uncertainty is 0
            % interets gives indication on how likely a target will stay at a grid
            obj = initProbability(obj);
            obj = initUncertainty(obj);
            obj = initInterests(obj);
            
            % initialise Receding Horizon Control parameters
            obj = initRHC(obj);
        end
        
        %% protected methods
        % ========================
        % get world centre
        % ========================
        function [obj] = getCentre(obj)
            obj.world_centre    = [obj.size_x/2, obj.size_y/2];
        end
        
        % ========================
        % get number of tiles
        % ========================
        function [obj] = getNumberOfTiles(obj)
            % +1 to the division to include the boundary
            obj.number_of_tiles_x = ceil(obj.size_x/obj.grid_size)+1;
            obj.number_of_tiles_y = ceil(obj.size_y/obj.grid_size)+1;
            obj.number_of_tiles   = obj.number_of_tiles_x*obj.number_of_tiles_y;
			
			obj.property_x_index = ceil(0.5*obj.buffer_size(1)/obj.grid_size)+1;
			obj.property_y_index = ceil(0.5*obj.buffer_size(2)/obj.grid_size)+1;
        end
        
        % ========================
        % get meshgrids
        % ========================
        function [obj] = getMeshGrids(obj)
            [obj.meshgrid_XX,obj.meshgrid_YY] = meshgrid(0:obj.grid_size:obj.size_x, 0:obj.grid_size:obj.size_y);
        end
        
        % ========================
        % get grid centres
        % ========================
        function [obj] = getGridCentres(obj)
            obj.grid_centre_x = obj.meshgrid_XX(1:end,1:end)+obj.grid_size/2;
            obj.grid_centre_y = obj.meshgrid_YY(1:end,1:end)+obj.grid_size/2;
        end
        
        % ========================
        % initialise probability
        % ========================
        function [obj] = initProbability(obj)
            obj.prob0               = 0.4;
            obj.prob_nom            = 0.5;
            obj.prob_tau            = 0.99939;
            obj.prob_detect         = 1.0;
            obj.prob_detect_tau     = 0.4;
            obj.prob_nodetect       = 0.0;
%             obj.prob_nodetect_tau   = 0.7;
            obj.prob_nodetect_tau   = 0.5;
            obj.prob_range          = 20; % This was 50m
            
            obj.prob_global         = ones(obj.number_of_tiles_y, obj.number_of_tiles_x).*obj.prob0;
        end
        
        % ========================
        % initialise uncertainty
        % ========================
        function [obj] = initUncertainty(obj)
            obj.uncert0             = 1;
            obj.uncert_nom          = 1;
            obj.uncert_tau          = 0.99;
            
            obj.uncert_global       = ones(obj.number_of_tiles_y, obj.number_of_tiles_x).*obj.uncert0;
        end
        
        % ========================
        % initialise interests
        % ========================
        function [obj] = initInterests(obj)
            obj.inter0              = 1;
            obj.inter_wild          = 0.1;
            obj.inter_nom           = 1;
            obj.inter_tau           = 0.9999;
            obj.inter_nogo_radius   = 20; % was 50m % The higher this number the more the agent will be in standby applying pressure at its maximum radius since it cannot reach the actual target
            obj.inter_alert_radius  = 50; % was 450m
            obj.inter_nogo          = 0;
            
            obj.inter_global        = ones(obj.number_of_tiles_y, obj.number_of_tiles_x).*obj.inter_wild;
			tmp_x = obj.property_x_index;
			tmp_y = obj.property_y_index;
            obj.inter_global(tmp_y:end-tmp_y,tmp_x:end-tmp_x) = obj.inter_global(tmp_y:end-tmp_y,tmp_x:end-tmp_x).*(obj.inter0/obj.inter_wild);
        end
        
        % ========================
        % initialise Receding Horizon Control parameters
        % ========================
        function [obj] = initRHC(obj)
            obj.RHC_steps           = 10;%22;
        end
        
        %% private methods
        % ========================
        % expand maps
        % ========================
        function [obj] = expandArrays(obj, iterations)
            obj.prob_global(:,:,2:iterations)   = zeros(obj.number_of_tiles_y, obj.number_of_tiles_x, iterations-1);
            obj.uncert_global(:,:,2:iterations) = zeros(obj.number_of_tiles_y, obj.number_of_tiles_x, iterations-1);
            obj.inter_global(:,:,2:iterations)  = zeros(obj.number_of_tiles_y, obj.number_of_tiles_x, iterations-1);
        end
        
        % ========================
        % update global interests
        % ========================
        function [obj] = updateGlobalInterests(obj, targets, cumulation)
            % loop through all the targets to get their interests map
            interests_global_temp               = zeros(obj.number_of_tiles_y, obj.number_of_tiles_x, size(targets,1));
            for k = 1:size(targets,1)
                interests_global_temp(:,:,k)    = targets(k).inter_map(:,:,cumulation);
            end
            % global interests is the minimum of all targets
            obj.inter_global(:,:,cumulation)    = min(interests_global_temp,[],3);
        end
        
        % ========================
        % update global uncertainty
        % ========================
        function [obj] = updateGlobalUncertainty(obj, agents, cumulation)
            % loop through all the agents to get their interests map
            uncertainty_global_temp             = zeros(obj.number_of_tiles_y, obj.number_of_tiles_x, size(agents,1));
            for i = 1:size(agents,1)
                uncertainty_global_temp(:,:,i)  = agents(i).uncert_map(:,:,cumulation);
            end
            % global uncertainty is the minimum of all agents
            obj.uncert_global(:,:,cumulation)   = min(uncertainty_global_temp,[],3);
        end
        
        % ========================
        % update global probability
        % ========================
        function [obj] = updateGlobalProbability(obj, agents, cumulation)
            % loop through all the agents to get their detection map
            intersection_detect_global = false(obj.number_of_tiles_y, obj.number_of_tiles_x);
            intersection_nodetect_global = false(obj.number_of_tiles_y, obj.number_of_tiles_x);
            for i = 1:size(agents,1)
                intersection_detect_global = intersection_detect_global|agents(i).detect_map;
                intersection_nodetect_global = intersection_nodetect_global|agents(i).nodetect_map;
            end
            % extract data to make code less cluttered
            data    = obj.prob_global(:,:,cumulation-1);
            data    = data.*obj.prob_tau + obj.prob_nom*(1-obj.prob_tau);
            data(intersection_nodetect_global) = data(intersection_nodetect_global).*obj.prob_nodetect_tau + obj.prob_nodetect*(1-obj.prob_nodetect_tau);
            data(intersection_detect_global) = data(intersection_detect_global).*obj.prob_detect_tau + obj.prob_detect*(1-obj.prob_detect_tau);
            obj.prob_global(:,:,cumulation) = data;
        end
    end
end
            