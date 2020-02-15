classdef VolumeDiscretization < matlab.mixin.Copyable
    % VolumeDiscretization - (Un)Occupancy map stored as dense 3D matrix indicating interior voxels
    %
    % We use three indexing schemas for the data: coordinates, offsets, and index
    %
    % # Coordinates represent the original, continuous, data space. We assume each point is represented by the tuple
    % (x,y,z). Coordinates are associated with voxel centers.
    %
    % # Offsets represent the zero-based, ordered, axis offsets associated with the quantized spatial coordinate. This is
    % also a three-tuple, however the ordering is affected by the property axisorder.
    %
    % # Index is a scalar index into the contiguous matrix V. This is the linear index associated with vect(V).
    %
    %
    properties
        V;      % Volume: dense matrix ordered according to axisorder
    end
    properties(SetAccess = protected)
        deltav; % voxel width (spacing)
        origin; % 3-vector identifying spatial location of V(1) voxel center in x,y,z coordinates
        grow;   % scalar - number of extra voxels to add when growing volume
    end
    
    properties(SetAccess = protected)
        axisorder; % Ordering of x,y,z axes in V. E.g. [2,3,1] specifies y-major, then z and finally x
    end
    properties(Dependent)
        axisorderinverse; % Ordering of x,y,z axes in V. E.g. [2,3,1] specifies y-major, then z and finally x
        medianpt;
    end
    methods
        function obj = VolumeDiscretization(V,varargin)
            % Constructor
            %
            %  vd = VolumeDiscretization(V);
            %  vd = VolumeDiscretization(sz); where sz is a 3-element vector of positive integers
            %  vd = VolumeDiscretization(...,origin); where sz is a 3-element vector of positive integers
            %  vd = VolumeDiscretization(...,deltav); deltav is a scalar voxel width
            %
            
            p = inputParser;
            default_order = [2,1,3];
            addOptional(p,'origin',zeros(3,1), @(x) validateattributes(x,{'numeric'},{'numel',3}));
            addOptional(p,'deltav',1, @(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
            addParameter(p,'order',default_order, @(x)isequal(sort(double(x(:))),[1;2;3]));
            addParameter(p,'grow',10, @(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
            
            parse(p,varargin{:});
            
            obj.axisorder = p.Results.order;
            obj.deltav = p.Results.deltav;
            obj.origin = p.Results.origin(:);     % Ensure column-width
            obj.grow = p.Results.grow;
            
            if nargin < 1
                V = repmat(100,1,3);
            end
            if  numel(size(V))~=3
                obj.V = zeros(V,'int8');    % no need for floats
            else
                obj.V = V;
            end
        end
        
        function [varargout] = isosurface(obj,varargin)
            % isosurface - Call build-in isosurface using X,Y,Z,V from VolumeDiscretization
            
            % origin is in coordinates with axis order (x,y,z). However, V is in axisorder. Plot everything in x,y,z
            % order.
            
            % isosurface of volume
            xyz = cell(1,3);
            % Re-order [x,y,z] rows in axisorder
            aoi = obj.axisorder;
            dims = size(obj.V);
            dims = dims(aoi)-1;
            for idx = 1:3
                xyz{idx} = obj.origin(idx):obj.deltav:obj.origin(idx)+obj.deltav*dims(idx);
            end
            
            varargout = cell(1,nargout);
            [varargout{:}] = isosurface(xyz{:},obj.V,varargin{:});            
        end
        
        function [varargout] = scatter3(obj,varargin)
            % scatter3
            Y = obj.index2coordinate(find(obj.V)); %#ok
            varargout = cell(1,nargout);
            [varargout{:}] = scatter3(Y(1,:),Y(2,:),Y(3,:),varargin{:});
        end
        
        function cover(obj,X)
            % cover - expand volume to cover input coordinates
            %
            %  X - 3xN matrix of coordinates
            %
            
            validateattributes(X,{'numeric'},{'nrows', 3});
            
            % Process lower then upper bound
            sgn = [-1,1];
            % Shift extra points to third dim
            X = permute(X,[1,3,2]);
            % Find required extension in two directions : res should be 3x2 matrix
            res = max((X - obj.bounds()).*sgn,[],3);
            % For each bound [lower,upper]
            for bidx = 1:numel(sgn)
                % Look for axes to extend
                dim = find(res(:,bidx)>0);
                if isempty(dim)
                    continue;
                end
                for ax = dim(1:numel(dim)).'
                    % Convert distance to samples and add grow term
                    ext = ceil(res(ax,bidx)/obj.deltav) + obj.grow;
                    extend(obj,ext*sgn(bidx),ax);
                end
            end
        end
        
        function rslt = get.axisorderinverse(obj)
            % Invert axis order
            [~,rslt] = sort(obj.axisorder);
        end
        
        function vd = subblock(obj,X0,Xend)
            % Return new VolumeDiscretization representing a subblock of the existing data.
            %
            % Input coordinates will be quantized (to the nearest voxel center). The data will not be sampled, only
            % shifted and truncated.
            %
            X = [X0(:),Xend(:)];
            if ~isempty(find(diff(X,1,2)<0,1))
                error('Inputs must bound volume');
            end
            offset = obj.coordinate2offset(X);
            new_origin = obj.offset2coordinate(offset(:,1));
            offset = offset +1;
            vd = VolumeDiscretization(obj.V(offset(1,1):offset(1,2),offset(2,1):offset(2,2),offset(3,1):offset(3,2)),...
                new_origin,obj.deltav,'order',obj.axisorder);
        end
        
        function idx = coordinate2index(obj,X)
            % coordinate2index - Discretize x,y,z coordinates to index
            %
            %  idx = obj.coordinate2index(X)
            %
            %  X : 3xN matrix of coordinates
            %  idx : 1xN vector of voxel indices
            %
            Y = coordinate2offset(obj,X);
            idx = offset2index(obj,Y);
        end
        
        function X = index2coordinate(obj,idx)
            % index2coordinate - Convert index to  x,y,z coordinates
            %
            %  X = obj.index2coordinate(idx)
            %
            %  idx : 1xN vector of voxel indices
            %  X : 3xN matrix of coordinates
            %
            Y = index2offset(obj,idx);
            X = offset2coordinate(obj,Y);
        end
        
        function Y = coordinate2offset(obj,X)
            % coordinate2offset - Discretize x,y,z coordinates to offsets
            %
            %  Y = obj.coordinate2offset(X)
            %
            %  X : 3xN matrix of coordinates
            %  Y : 3xN vector of integer-valued offsets
            %
            validateattributes(X,{'numeric'},{'nrows', 3});
            
            % Convert to integer-values
            Y = round((X - obj.origin)/obj.deltav);
            Y = Y(obj.axisorderinverse,:);
            Y = max(min(Y,size(obj.V).'-1),0);
        end
        
        function idx = offset2index(obj,Y)
            % offset2index - Convert axis offsets to voxel indices
            %
            %  idx = obj.offset2index(X)
            %
            %  Y   : 3xN vector of integer-valued offsets
            %  idx : 1xN vector of voxel indices
            %
            validateattributes(Y,{'numeric'},{'integer','nrows', 3});
            
            % Offset for subsequent dimensions
            dims = size(obj.V);
            idx = 1+Y(1,:) + cumprod(dims(1:2))*Y(2:end,:);
        end
        
        function Y = index2offset(obj,idx)
            % index2offset - Convert voxel indices to axis offsets
            %
            %  Y = obj.index2offset(idx)
            %
            %  idx : 1xN vector of voxel indices
            %  Y   : 3xN vector of integer-valued offsets
            %
            dims = size(obj.V);
            tmp = cumprod(dims(1:2));
            Y = zeros(3,numel(idx));
            idx = idx(:).' -1;
            for itt = 3:-1:2
                Y(itt,:) = floor(idx/tmp(itt-1));
                idx = idx - tmp(itt-1)*Y(itt,:);
            end
            Y(1,:) = idx;
        end
        
        function X = offset2coordinate(obj,Y)
            % offset2coordinate - Convert axis offsets to coordinates of pixel centers
            %
            %  X = obj.index2offset(Y)
            %
            %  Y : 3xN vector of integer-valued offsets
            %  X : 3xN matrix of coordinates
            validateattributes(Y,{'numeric'},{'integer','nrows', 3});
            
            X = Y(obj.axisorder,:)*obj.deltav + obj.origin;
        end
        
        function bnds = bounds(obj)
            % Bounds: return voxel centers at origin and extrema
            %
            %  result : 3x2 matrix. First column is the origin.
            dims = size(obj.V).';
            bnds = [obj.origin, obj.origin+(dims(obj.axisorder)-1)*obj.deltav];
        end
        
        function extend(obj,N,ax_xyz)
            % Extend the volume N voxels along specified axis
            %
            % For N < 0 , we assume the origin will change and while adding |N| voxels
            %
            
            ax = obj.axisorder(ax_xyz);
            
            % Construct fill
            tmp = size(obj.V);
            tmp(ax) = abs(N);
            fill = zeros(tmp);
            
            if N > 0
                % Positive extension
                %  zero pad
                %  no change to origin
                obj.V = cat(ax,obj.V,fill);
            else
                % Negative extension
                %  concatenate zeros
                %  origin will be updated
                obj.V = cat(ax,fill,obj.V);
                obj.origin(ax_xyz) = obj.origin(ax_xyz) + N*obj.deltav;
            end
        end
        
        function permute(obj,p)
            % Permute - reorder dimensions of p
            obj.V = permute(obj.V,p);
            obj.deltav = obj.deltav(p);
            obj.origin = obj.origin(p);
            obj.axisorder = obj.axisorder(p);
        end
        
        function rslt = and(obj,varargin)
            % BITAND - Only consider V
            rslt = copy(obj);
            rslt.V = 1*and(obj.V,varargin{:});
        end
        
        function val = get.medianpt(obj)
            val = size(obj.V).';
            val = val(obj.axisorderinverse);
            val = floor(0.5*val) + obj.origin;
        end
        
    end
    methods(Static)
        function num = axis(char)
            % Convert x,y,z to 1,2,3
            num = upper(char) - 'W';
        end
    end
end

