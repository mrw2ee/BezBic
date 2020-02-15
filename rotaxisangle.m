function [R,varargout] = rotaxisangle(u,theta,flag)
% ROTAXISANGLE - Rotation matrix from axis-angle
%
% R = rotaxisangle(u,theta) for 3-vector u, and angle theta in radians.
% R = rotaxisangle(u,x,'uniform') supply uniformly distributed x [0,1].
%	In this case theta will be derived from the cumulative distribution function (theta - sin theta)/pi.
%   This ensures the 3D rotations are uniformly distributed.
%
% [R,theta] = rotutheta(u,x,'uniform') also returns the rotation angle theta.
%
% Construct the rotation matrix R such that
%  y = Rx
%
% Note, u is NOT NORMALIZED internally. As such, the identify matrix is
% returned for 
%  I = rotutheta([1 1 1],0)
%
% M.Walker 5/17/2017

if nargin > 2
    if isequal(flag,'uniform') && theta>=0 && theta<=1
        theta = fzero(@(x)pi*theta-(x-sin(x)),pi/2);
    else
        error('Third argument must be ''uniform'' and 0 <= theta <=1');
    end
end
u = reshape(u,[],1);
sinth = sin(theta);
costh = cos(theta);

R = [costh+u(1)^2*(1-costh), u(1)*u(2)*(1-costh)-u(3)*sinth, u(1)*u(3)*(1-costh)+u(2)*sinth;...
    u(2)*u(1)*(1-costh)+u(3)*sinth, costh+u(2)^2*(1-costh), u(2)*u(3)*(1-costh)-u(1)*sinth;...
    u(3)*u(1)*(1-costh)-u(2)*sinth, u(3)*u(2)*(1-costh)+u(1)*sinth, costh + u(3)^2*(1-costh)];

if nargout > 1
    varargout{1} = theta;
end

