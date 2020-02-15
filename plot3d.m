function varargout = plot3d(fname,varargin)
% PLOT3D - Setup figure for 3D model with lighting

%% Setup figure
fhdl = gcf;
% Camera toolbar does not clean up after itself
%clf(fhdl,'reset');
% Rotate3D only allows 180 degrees in elevation changes. Cameratoolbar does
% not suffer this limitation.
cameratoolbar(fhdl,'SetMode','orbit');
cameratoolbar(fhdl,'Show');
% Map zoom to scroll wheel
set(fhdl,'WindowScrollWheelFcn',@(~,evd)camzoom(1.15^-evd.VerticalScrollCount));

%% Generate plot
varargout = cell(1,nargout);
[varargout{:}] = fname(varargin{:});


%% Additional setup 
% Light model
lightwithcam;
lighting gouraud
material dull

% Now that we have axes, this will stick
cameratoolbar(fhdl,'SetCoordSys','none');

% Turn Orbit mode 'off' when user pressed any key on keyboard. 
%  e.g. allowing the user to select a new point on the surface.
set(gcf,'WindowKeyPressFcn',@(src,evt)cameratoolbar(src,'SetMode','nomode'));