classdef VolumeAccumulator < handle
    % VolumeAccumulator - Track "visited" space of catheter
    %
    % Maintains a map of the interior volume (property .vd) and adds a method to process new position data
    %
    
    properties
        vd;                     % VolumeDiscretization
        cath_discretization;    % catheter length discretization. Will be multiplied with orientation
        cath_radius;            % catheter radius
    end
    
    methods
        function obj = VolumeAccumulator(varargin)
            % VolumeAccumulator Constructor
            %
            % VolumeAccumulator(vd) - Specify VolumeDiscretization
            % VolumeAccumulator(...,'cath_discretization',d) - Vector specifying axial offset
            % VolumeAccumulator(...,'cath_radius',r) - scalar radius of catheter
            % 
            p = inputParser;
            
            addOptional(p,'vd',VolumeDiscretization(), @(x) validateattributes(x,{'VolumeDiscretization'},{'scalar'}));
            addParameter(p,'cath_discretization',(-10:0), @(x) validateattributes(x,{'Numeric'},{'vector'}));
            addParameter(p,'cath_radius',1, @(x) validateattributes(x,{'Numeric'},{'vector'}));
            
            parse(p,varargin{:});
            
            obj.vd = p.Results.vd;
            obj.cath_discretization = p.Results.cath_discretization;
            obj.cath_radius = p.Results.cath_radius;
        end
        
        function accumulate(obj,tip,ori)
            % accumulate - Flag voxels in volumediscretization
            %
            % tip - 3xN matrix of tip locations
            % ori - 3xN matrix of unit vectors specifying tip orientation at each sample
            if nargin == 2
                validateattributes(tip,{'numeric'},{'nrows', 6});
                ori = tip(4:6,:);
                tip(4:6,:) = [];
            else
                validateattributes(tip,{'numeric'},{'nrows', 3});
                validateattributes(ori,{'numeric'},{'nrows', 3});
            end
                
            if size(tip,2)~= size(ori,2)
                error('Inputs must have same number of samples');
            end
            
            % Additional offsets due to catheter tip - axial only
            %   Third dimension will represent discretization points
            cath_axial = reshape(kron(obj.cath_discretization,ori),[size(ori),numel(obj.cath_discretization)]);
            cath_coords = tip + cath_axial;
            %cath_coords = tip;
            
            % Should we extend the volume?
            pad = 5*obj.vd.deltav; % Minimum pad for extra points
            % Endpoint are sufficient
            tmp = [cath_coords(:,:,1),cath_coords(:,:,end)];
            tmp = [min(tmp,[],2)-pad, max(tmp,[],2)+pad];
            obj.vd.cover(tmp);
            
            % Convert all coordinates to scalar indices
            idx = unique(obj.vd.coordinate2index(reshape(cath_coords,3,[])));
            %%{
            % Offsets due to finite catheter radius
            dims = size(obj.vd.V);
            rad_bins = ceil(obj.cath_radius/obj.vd.deltav);
            idx = idx + [0;kron([-rad_bins:-1,1:rad_bins],[1,cumprod(dims(1:2))]).'];
            %}
            % Turn on all touched voxels
            obj.vd.V(idx) = 1;
        end
    end
end