function R = pointselection(V,varargin)
% POINTSELECTION - Indicate surface voxels and weighted inverse distance from center
%
% R = pointselection(V)
%   V - Dense matrix indicating interior voxels
%   R - Dense matrix of weights approximating inverse distance of surface voxels to the center.
%
% pointselection(...,'rgrow',r)
%   r : positive scalar indicating the number of times to 'grow' the surface
%
% pointselection(...,'eps',eps)
%   eps: positive scalar specifying the minimum number of exterior voxels in a neighborhood to indicate boundary
%
% We identify surface voxels near the center of the volume. Iteratively we expand the search in the neighborhood of
% known surface voxels. The iteration index approximates the distance, along the surface, from new voxels to the
% starting point. For R(n)>0, the iteration index is given by 1+max(R)-R(n). In this way, R(n)\in[0,M], approximates a
% weighted inverse distance from center.
%
% M.Walker 10/17/2019

persistent p
if isempty(p)
    default_rgrow = 2;
    default_eps = 9;
    p = inputParser;
    addRequired(p,'V',@(x)validateattributes(V,{'numeric'},...
        {'3d','nonempty'}))
    addParameter(p,'rgrow',default_rgrow,@(x)validateattributes(x,...
        {'numeric'},{'scalar','integer'}))
    addParameter(p,'eps',default_eps,@(x)validateattributes(x,...
        {'numeric'},{'scalar','nonnegative'}))
end
parse(p,V,varargin{:})
rgrow = p.Results.rgrow;

% Verify input is of sufficient size
N = size(V);
if ~isempty(find(N < 2*rgrow+1,1))
    error('Input volume too small');
end
% Number of times to iteratively blur
% We assume boundary voxels (of max dim) are ambiguous. So, stop before we get to that point.
ngrow = floor((0.5*(max(N)-1)-1)/rgrow);

vtype = class(V);
% Laplacian kernel
W = -ones([3,3,3],'like',V);
W(2,2,2) = 3^3-1;
% Determine boundary voxels
V = cast(convn(V,W,'same')>=p.Results.eps,vtype);

% Starting blur indicator
R = zeros(N,vtype);
idx = arrayfun(@(n)ceil(0.5*n)-rgrow:ceil(0.5*n)+rgrow,N,'UniformOutput',false);
R(idx{:}) = V(idx{:});

for itt = 2:ngrow
    tmp = convn(R,ones(repmat(2*rgrow+1,1,3),vtype),'same');
    R = R + cast(tmp & V,vtype);
end
