% BezBic Demo
%
% Visualize online surface approximation on procedure data
%
% This is a challening problem highlighting the limits of our approach. Other demos highlight the benefits.
%
% For procedure data the occupancy map contains many missing samples where interior voxels are not traversed by the
% catheter. When the occupancy map is inaccurate near the tip of the catheter, local surface approximation is
% challenging.
%
% This script is basically an animation. The user can rotate the anatomic map by clicking on the plot and zooming in
% with the scroll wheel. A good way to pause the animation is placing a breakpoint in the for loop.
%
% Author: Michael R. Walker II <mwalkerii@wustl.edu>
%
% This file is part of the BezBic library. See 'LICENSE' for copyright and citation info.

%% Setup
% Algorithm and simulation parameters
clear;

%%% Simulation parameters
% Procedures typically begin with a mapping phase. During the mapping phase the occupancy map is inaccurate,
% particularly near the tip of the catheter. The current local surface approximation algorithms focus on more reliable
% occupancy maps. So, we start late in the procedure.
n0 = 5e4;
% Specify number of position samples to buffer between updates (positions sampled at 7.5 Hz)
Nstep = 7;

%%% Algorithm parameters
eps = 9;            % Threshold for boundary voxel (minimum number of exterior voxels in 3-cubic volume)
lambda = 1e-2;      % Regularization on Bezier surface control points.

% Optimization parameters specifying stopping criteria, etc.
%  optsA is used by fitsurface (and fitsurfacefixedorder)
%  optsuv is used by fminnewton
optsA = struct('MaxIter',10,'TolX',1e-3,'TolFun',1e-2,'lambda',lambda);
optsuv = struct('MaxIter',10,'TolX',1e-3,'TolFun',1e-2,'MaxFunEvals',20,'sigma',1e-1,'beta',0.1);

%%% Volume extraction:
% Working with the full occupancy map is computationally expensive, particularly when performing 3D convolutions. For
% many steps we extract a small cubic volume and perform computations locally. We specify the size of the volumes as
% (2*vol+1)^3

% Initial point correction: Often the catheter tip is not in contact with the surface. The tip position is not a
% suitable starting point for surface fitting. Search for the nearest boundary voxel using an extracted volume.
vol_correction = 5;
% Point selection: Selecting points approximating a single surface require multiple convolutions.
vol_pointselect = vol_correction+1;

%% Load data
% Load procedure data and initialize occupancy map

% Positions: 6 vector indicating catheter tip position and orientation. Position is in mm
load('ProcData03','positions');

% Initialize occupancy map
va = VolumeAccumulator;
va.accumulate(positions(:,1:n0-1));
n = n0;

%% Initialize plots

% Find local point on surface
pt = positions(1:3,n-1);
ptx = pointcorrect(va.vd,pt,eps,vol_correction);
% Get a point cloud approximating surface (using pointselection)
Y = pointcloud(va.vd,ptx,eps,vol_pointselect);
% Fit local surface
[Ahat,uvhat,sigma2hat,exitflag,rslt] = fitsurface(Y,1,optsA,optsuv);
X = meshsurf(Ahat);

figure(1)
clf;
% Plot occupancy map
fv = va.vd.isosurface();
hmap = plot3d(@patch,fv,'EdgeColor','none','FaceColor',[0.5 .2 .2],'FaceAlpha',0.5);
material dull
hold on

% Plot mesh of parametric surface
xdim = sqrt(size(X,2));
dat = X.';
hsurf = surf(reshape(dat(:,1),xdim,xdim),reshape(dat(:,2),xdim,xdim),reshape(dat(:,3),xdim,xdim),'FaceColor',[0.4 0.4 1],'FaceAlpha',0.5);

% Current point selection
hpoints = scatter3(Y(1,:),Y(2,:),Y(3,:),'g.','SizeData',200,'LineWidth',2);
hpt = scatter3(pt(1,:),pt(2,:),pt(3,:),'o','SizeData',100,'LineWidth',2);
hptx = scatter3(ptx(1,:),ptx(2,:),ptx(3,:),'mx','SizeData',100,'LineWidth',2);

set(gca,'CameraTarget',ptx);
legend('Occupancy Map','BezBic Surface','Point Cloud','Tip Position','Starting Point');

%%
K = floor((size(positions,2)-n)/Nstep);
for k = 1:K
    
    % Feed occupancy map
    va.accumulate(positions(:,n:n+Nstep-1));
    n = n+Nstep;
    
    % Find local point on surface
    pt = positions(1:3,n-1);
    ptx = pointcorrect(va.vd,pt,eps,vol_correction);
    
    % Algorithm  ------------------------------------------------------------------------------------------------------
    % Get a point cloud approximating surface (using pointselection)
    Y = pointcloud(va.vd,ptx,eps,vol_pointselect);
    % Fit local surface (when the catheter is near a surface)
    if ~isempty(Y)
        [Ahat,uvhat,sigma2hat,exitflag,rslt] = fitsurface(Y,1,optsA,optsuv);
    end
    %% End of algorithm -------------------------------------------------------------------------------------------------
    
    % Generate local isosurface of occupancy map. This is nonparametric and may be messy
    vol = 20;
    v0 = va.vd.subblock(pt-vol,pt+vol);
    fv = v0.isosurface();   % Note, this is computationally expensive
    % Generate mesh of parametric surface for visualization
    X = meshsurf(Ahat);
    xdim = sqrt(size(X,2));
    dat = X.';
    
    % Update plots
    set(hmap,'Faces',fv.faces,'Vertices',fv.vertices);
    set(hsurf,'XData',reshape(dat(:,1),xdim,xdim),'YData',reshape(dat(:,2),xdim,xdim),'ZData',reshape(dat(:,3),xdim,xdim));
    set(hpoints,'XData',Y(1,:),'YData',Y(2,:),'ZData',Y(3,:));
    set(hpt,'XData',pt(1),'YData',pt(2),'ZData',pt(3));
    set(hptx,'XData',ptx(1),'YData',ptx(2),'ZData',ptx(3));
    drawnow;
    
    % <-------- Good spot for a breakpoint
    if (mod(k,20)==0)
        % Smoothly update camera target as catheter traverses the anatomy
        updatecam(gca,ptx);
    end
end


%% Local functions
function ptx = pointcorrect(vd,pt,eps,vol)
% pointcorrect - Move point to nearest boundary point
%
% ptx = pointcorrect(vd,pt,eps,vol)
%  vd - VolumeDiscretization object
%  pt - 3x1 vector of initial point
%  eps - Threshold for identifying boundary voxels with Laplacian kernel
%  vol - Limit the search volume about pt
%  ptx - Corrected point (on surface)
%
% Note, if no boundary points in extracted volume, returns empty

% Extract smaller volume
v0 = vd.subblock(pt-vol,pt+vol);
V = v0.V;

% Laplacian kernel
W = -ones([3,3,3],'like',V);
W(2,2,2) = 3^3-1;
% Determine boundary voxels
idx = find(convn(V,W,'same')>=eps);
if isempty(idx)
    ptx = zeros(3,0,'like',pt);
else
    % Convert to coordinates
    Y = v0.index2coordinate(idx);
    % Return closest nearest
    [~,idx] = min(sum((Y-pt).^2));
    ptx = Y(:,idx);
end
end



function Y = pointcloud(vd,pt,eps,vol)
% pointcloud - Identify point approximating single surface
%
% Y = pointcloud(vd,pt,eps,vol)
%  vd - VolumeDiscretization object
%  pt - 3x1 vector of initial point
%  eps - Threshold for identifying boundary voxels with Laplacian kernel
%  vol - Limit the search volume about pt
%  Y - 3xN point cloud

% Extract smaller volume
v0 = vd.subblock(pt-vol,pt+vol);

R = pointselection(v0.V,'eps',eps);
% Convert output
idx = find(R>0);
Y = v0.index2coordinate(idx);
end

function [X,varargout] = meshsurf(A,n)
% MESHSURF - Mesh Bezier surface with control points A
%  X = meshsurf(A)
%  X = meshsurf(A,n) - Optionally specify the number of patches in each surface dim
%  [X,UV] = meshsurf(A) - Optionally return the uv coordinates: 2 x (n+1)^2 matrix
%
if nargin < 2
    n = 10;
end
[U,V] = meshgrid(linspace(0,1,n+1),linspace(0,1,n+1));
X = bezsurf.uv2x(A,[U(:),V(:)].');
if nargout>1
    varargout = {[U(:),V(:)].'};
end
end

function updatecam(htarget,pt)
% Smoothly update camera target to new location
pt0 = camtarget(htarget);
d = pt.' - pt0;
N = 100; % Number of update steps
period = 16e-3; % 60 Hz
fctn = @(s,~,h)camtarget(h,pt0+d*0.5*(1-cos(pi*s.TasksExecuted/N)));
tmr = timer('TimerFcn',{fctn,htarget},'Period',period,'TasksToExecute',N,'ExecutionMode','fixedRate');
start(tmr);
end