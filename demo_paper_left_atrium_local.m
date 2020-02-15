%% Script visualizing procedure data and local surface fitting
%
% Author: Michael R. Walker II <mwalkerii@wustl.edu>
%
% This file is part of the BezBic library. See 'LICENSE' for copyright and citation info.


%% Setup

clear;

% Load data and construct occupancy map
load('ProcData03','positions','timestamps','surfaces');
va = VolumeAccumulator;
% Read all positions at once
va.accumulate(positions);

% for saving figures...
output_dir = "..\doc\";
file_prefix = "plt_leftatrium_";
file_prefix = "";
plt = struct('dims',[3.4,2.8]);
font_sz = 8;

figure('Units','inches','DefaultAxesFontSize',font_sz);
clf;
hlspv = [];
hlipv = [];
hlaa = [];

% Plot TRI surface - interpolated
hsurf = plot3d(@trisurf,surfaces{1},'EdgeColor','none','FaceColor',[0.4 0.4 1],'FaceAlpha',0.5);
hsurf.PickableParts = 'none';
set(gca,'FontSize',font_sz);
hold on

% Plot occupancy map
fv = va.vd.isosurface();
hpos = patch(fv,'EdgeColor','none','FaceColor',[0.5 .2 .2],'FaceAlpha',0.5);
material dull

% Current point selection
h = scatter3(0,0,0,'.g','SizeData',100,'LineWidth',2);
set(hpos,'ButtonDownFcn',{@update_xyz,va.vd,h});


%% Visualize point selection

hl = legend([hsurf,hpos,h],{'mesh','$V$','$\mathcal{X}$'},'Interpreter','Latex','Location','SouthWest');

% Anatomic landmarks
if ~isempty(hlspv)
    delete(hlspv);
    delete(hlipv);
    delete(hlaa);
end
lspv = [105 -121 107];
lipv = [114 -164 140];
laa = [108 -114 118];
hlspv = text(lspv(1),lspv(2),lspv(3),'LSPV','HorizontalAlignment','left','FontSize',font_sz);
hlipv = text(lipv(1),lipv(2),lipv(3),'LIPV','HorizontalAlignment','right','FontSize',font_sz);
hlaa = text(laa(1),laa(2),laa(3),'LAA','HorizontalAlignment','left','FontSize',font_sz);

parms = struct('pt0',[85.5226 -128.8076 115.5050],...
    'CameraTarget',[96.5411 -126.7368 112.2764],...
    'CameraUpvector', [0.5476 -0.6287 -0.6200],...
    'CameraPosition', [-343.6780 335.2564 -641.7989],...
    'CameraViewAngle',2.0524);

% Select points
src = hpos;
evt = struct('IntersectionPoint',parms.pt0);
update_xyz(src,evt,va.vd,h);    % Update plot

% Update plot
set(gca,'CameraTarget',parms.CameraTarget,'CameraPosition',parms.CameraPosition,...
    'CameraUpvector',parms.CameraUpvector,'CameraViewAngle',parms.CameraViewAngle);

plt_name = "data";
if strlength(file_prefix)
    % Squeeze the plot and save *.eps
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2),plt.dims]);
    % Expand axis to figure (no pad for labels)
    %set(gca,'Position',[0 0.05 1 0.95])
    set(hl,'Position',[0,0,hl.Position(3:end)]);
    axis image
    axis off
    print('-r300',[output_dir+file_prefix+plt_name],'-depsc');
end

%% Local surface approximation
%----------------------------------------------------------------------------------------------------------------------
%  This is the full algorithm without plotting

% Get point cloud
pt = parms.pt0; 
Y = pointcloud(va.vd,pt);

% Fit local surface
lambda = 1e-4;
optsA = struct('MaxIter',10,'TolX',1e-5,'TolFun',1e-5,'lambda',lambda);
optsuv = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',20,'sigma',1e-1,'beta',0.1);
[Ahat,uvhat,sigma2hat,exitflag,rslt] = fitsurface(Y,1,optsA,optsuv);

%----------------------------------------------------------------------------------------------------------------------
% Get footpoints
Xhat = bezsurf.uv2x(Ahat,uvhat);

fprintf('sigma2hat: %f\n',sigma2hat)


%% Find closest points on mesh
% Score results, and extract relevant triangles in surface mesh

% Compute distances
spts = surfaces{1}.Points;
d = zeros(1,size(Y,2));
pidx = d;
for n = 1:numel(d)
    [d(n),pidx(n)] = min(sum((spts - Y(:,n).').^2,2));
end
d = sqrt(d);
fprintf('Median distance to nonparametric: %f\n',median(d))

% Construct abbreviated mesh
s2 = surfaces{1};
valid_pt = false(size(s2.ConnectivityList,1),1);
parfor n = 1:numel(valid_pt)
    valid_pt(n) = ~isempty(intersect(s2.ConnectivityList(n,:),pidx));
end
s2 = triangulation(s2.ConnectivityList(valid_pt,:),s2.Points);


%% Visualize local approximation

parms = struct('pt0',[85.5226 -128.8076 115.5050],...
    'CameraTarget',[85.5226 -128.8076 115.5050],...
    'CameraUpvector', [0.4829 0.7091 -0.5139],...
    'CameraPosition', [941.1453 -431.2831 502.1522],...
    'CameraViewAngle',0.6400);

figure('Units','inches','DefaultAxesFontSize',font_sz);
X = meshsurf(Ahat);
n = sqrt(size(X,2));
dat = X.';

% Plot section of nonparametric map
hsurfLocal = plot3d(@trisurf,s2,'EdgeColor','none','FaceColor',hsurf.FaceColor,'FaceAlpha',hsurf.FaceAlpha);
set(gca,'FontSize',font_sz);
hold on
% Plot mesh of parametric surface
s1 = surf(reshape(dat(:,1),n,n),reshape(dat(:,2),n,n),reshape(dat(:,3),n,n),'FaceColor',hpos.FaceColor,'FaceAlpha',hpos.FaceAlpha);

% Plot the selected points used to fit the surface
%  Additionally, color the markers by the distance to the surface approximation
dat = Y.';
d = sqrt(sum((Y-Xhat).^2));
fprintf('Median distance to parametric: %f\n',median(d))
scat1 = scatter3(dat(:,1),dat(:,2),dat(:,3),50,d.','x','LineWidth',1.5);
set(gca,'clim',[min(d),max(d)])

set(gca,'CameraTarget',parms.pt0,'CameraPosition',parms.CameraPosition,...
    'CameraUpvector',parms.CameraUpvector,'CameraViewAngle',parms.CameraViewAngle);
c = colorbar('eastoutside');
c.Label.String = "distance (mm)";
legend([hsurfLocal,s1,scat1],'mesh','$\hat{\mathcal{S}}$','$\mathcal{X}$','Interpreter','Latex','Location','SouthEast');
    
plt_name = "mesh";
if strlength(file_prefix)
    % Squeeze the plot and save *.eps
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2),plt.dims]);
    axis image
    axis off    
    print('-r300',[output_dir+file_prefix+plt_name],'-depsc');
end
%% Local functions -----------------------------------------------------------------------------------------------------

function update_xyz(src,evt,vd,h)
    % User-selected point
    pt = evt.IntersectionPoint;
    h.UserData = pt.';
    camtarget(src.Parent,pt);
    Y = pointcloud(vd,pt);
    set(h,'XData',Y(1,:),'YData',Y(2,:),'ZData',Y(3,:));
    % Turn orbit mode back on
    cameratoolbar(src.Parent.Parent,'SetMode','orbit');
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

if nargin < 3
    eps = 9;
    vol = 7;
end
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
