%% Script demonstrating point selection from position data
%
% This is an interactive demo allowing the user to update the local collection of points for surface fitting.
%
% Press spacebar to enable "target selection," then select a new target on the occupancy map.
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
 
% Initialize plot
figure(1);
clf;
hlspv = [];
hlipv = [];
hlaa = [];

% Plot TRI surface - interpolated
hsurf = plot3d(@trisurf,surfaces{1},'EdgeColor','none','FaceColor',[0.4 0.4 1],'FaceAlpha',0.5);
hsurf.PickableParts = 'none';
hold on

% Plot occupancy map
fv = va.vd.isosurface();
hpos = patch(fv,'EdgeColor','none','FaceColor',[0.5 .2 .2],'FaceAlpha',0.5);
material dull

% Current point selection
h = scatter3(0,0,0,'o','SizeData',100,'LineWidth',2);
set(hpos,'ButtonDownFcn',{@update_xyz,va.vd,h});


%% Label plot and initialize view

legend([hpos,hsurf,h],{'Visited Volume','Interpolated Surface','Point Cloud'},'Location','SouthWest')

% Anatomic landmarks
if ~isempty(hlspv)
    delete(hlspv);
    delete(hlipv);
    delete(hlaa);
end
lspv = [105 -121 107];
lipv = [114 -164 140];
laa = [108 -111 118];
hlspv = text(lspv(1),lspv(2),lspv(3),'LSPV','HorizontalAlignment','left');
hlipv = text(lipv(1),lipv(2),lipv(3),'LIPV','HorizontalAlignment','right');
hlaa = text(laa(1),laa(2),laa(3),'LAA','HorizontalAlignment','left');

% Initial point selection & view
pt0 = [85.5225906372070,-128.807556152344,115.504997253418];
src = hpos;
evt = struct('IntersectionPoint',pt0);
update_xyz(src,evt,va.vd,h);
parms = struct('pt0',[85.5226 -128.8076 115.5050],...
    'CameraTarget',[96.5411 -126.7368 112.2764],...
    'CameraUpvector', [0.5476 -0.6287 -0.6200],...
    'CameraPosition', [-343.6780 335.2564 -641.7989],...
    'CameraViewAngle',2.0524);
src = hpos;
evt = struct('IntersectionPoint',parms.pt0);
update_xyz(src,evt,va.vd,h);
set(gca,'CameraPosition',parms.CameraPosition,'CameraUpvector',parms.CameraUpvector,...
    'CameraViewAngle',parms.CameraViewAngle);

%% Local function - Point selection and graphics updates
function update_xyz(src,evt,vd,h)
    % User-selected point
    pt = evt.IntersectionPoint;
    h.UserData = pt.';
    
    % Smoothly update camera target to new location
    htarget = src.Parent;
    pt0 = camtarget(htarget);
    d = pt - pt0;
    N = 100; % Number of update steps
    period = 16e-3; % 60 Hz
    fctn = @(s,~,h)camtarget(h,pt0+d*0.5*(1-cos(pi*s.TasksExecuted/N)));
    tmr = timer('TimerFcn',{fctn,htarget},'Period',period,'TasksToExecute',N,'ExecutionMode','fixedRate');
    start(tmr);
        
    % Generate point cloud
    vol = 9;
    v0 = vd.subblock(pt-vol,pt+vol);
    % Indicate target point
    %v0.V(vol,vol,vol) = 1;
    R = pointselection(v0.V);
    idx = find(R>0);
    wts = R(idx);
    Y = v0.index2coordinate(idx);
    
    % Update scatter plot. Assign color by weights
    Cmap = h.Parent.Colormap;
    c = Cmap(ceil(double(wts)/double(max(wts))*size(Cmap,1)),:);
    set(h,'XData',Y(1,:),'YData',Y(2,:),'ZData',Y(3,:),'Cdata',c);
    
    % Turn orbit mode back on
    cameratoolbar(src.Parent.Parent,'SetMode','orbit');
end




