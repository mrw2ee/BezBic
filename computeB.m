function B = computeB(uv,nu,nv)
% COMPUTEB - Compute B matrix from product of Bernstein basis vectors
%
% Syntax
%
% B = computeB(uv,nu,nv)
%  uv - 2 x nx matrix
%  nu,nv - basis order (positive integers)
%  B - (nu+1)(nv+1) x nx

% Compute Bernstein basis vectors
% Sample the basis vectors at (u,v) points
if nu==nv
    tmp = bernsteinbasis(uv(:),nu);
    Bu = tmp(1:2:end,:).';
    Bv = tmp(2:2:end,:).';
else
    Bu = bernsteinbasis(uv(1,:).',nu).';
    Bv = bernsteinbasis(uv(2,:).',nv).';
end

% Form B
nx = size(uv,2);
B = reshape(permute(Bu,[1 3 2]).*permute(Bv,[3 1 2]),(nu+1)*(nv+1),nx);
