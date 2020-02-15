function [X] = uv2x(A,uv)
% uv2x - Evaluate the Bezier surface at querry points
%
% Format
%
% X = uv2x(A,uv)
% 
% A - control points: nu x nv x 3, or (nu^2) x 3 where
% uv - 2 x p query points
%
%
% X - 3 x p points on surface
%
% M.Walker 4/10/2019

validateattributes(uv,{'numeric'},{'nrows',2},2)

% Determine the number of basis vectors in U and V
Asz = size(A);
if Asz(end) ~= 3
    error('Last dim of A must be 3');
end
if numel(Asz)==3
    nu = Asz(1);
    nv = Asz(2);
else
    nu = sqrt(Asz(1));
    nv = nu;
    if mod(nu,1)~=nu
        error('Could not determine nu,nv');
    end
end

p = size(uv,2);

% Sample the basis vectors at (u,v) points
if nu==nv
    tmp = bernsteinbasis(uv(:),nu-1);
    Bu = tmp(1:2:end,:);
    Bv = tmp(2:2:end,:);
else
    Bu = bernsteinbasis(uv(1,:).',nu-1);
    Bv = bernsteinbasis(uv(2,:).',nv-1);
end

BuA = reshape(Bu*reshape(A,nu,3*nv),p,nv,3);
X = reshape(sum(BuA.*Bv,2),p,3).';
