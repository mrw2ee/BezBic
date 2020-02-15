function lightwithcam(varargin)
%% LIGHTWITHCAM - Move light with camera
%
% MATLAB's default lighting control requires independent manipulation of
% the camera and light sources. In contrast, most 3D modeling software
% moves the light source and camera together as a unit. This approach
% highlights structural subtleties of the model.
%
% lightwithcam - Creates light in current axes
% lightwithcam(ax) 
% lightwithcam(...) - Additional inputs passed to 'light'
%
% For a smooth surface with no specular reflection call:
%
% lighting gouraud;
% material dull;
%
%
% M.Walker 7/18/2016

if nargin > 0 && isa(varargin{1},'matlab.graphics.axis.Axes')
    ahdl = varargin{1};
    varargin(1) = [];
else
    ahdl = gca;
end

% Light the model: single light, colocated with camera
lhdl = light(ahdl,varargin{:});
lhdl = camlight(lhdl,'headlight');
bringLightToo = @(src,evnt)set(lhdl,'Position',get(ahdl,'CameraPosition'));
% Note, by default, camera moves reset lighting. So, we just follow up the
% lighting reset.
all_listeners = cell(1,4);
all_listeners{1} = addlistener(lhdl,'Position','PostSet',bringLightToo);
all_listeners{2} = addlistener(ahdl,'View','PostSet',bringLightToo);
all_listeners{3} = addlistener(ahdl,'CameraPosition','PostSet',bringLightToo);
all_listeners{4} = addlistener(ahdl,'CameraUpVector','PostSet',bringLightToo);
% Clean up after ourselves when the light is deleted
% For example, this will happen if the patches are deleted.
lhdl.DeleteFcn = @(~,~)cellfun(@delete,all_listeners);