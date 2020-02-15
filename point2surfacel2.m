function [c,grad,H] = point2surfacel2(uvk,x,A)
% point2surfacel2 - Objective function: l2 distance from point to Bezier surface coordinates
%
% [c,grad,H] = point2surfacel2(uvk,x,A)
%  uvk : 2-vector of surface coordinates
%  x : 3x1 vector 
%  A : nu+1 x nv+1 x 3 matrix of control points
%
%  c : Scalar (positive) cost
%  grad : 2x1 gradient
%  H : 2x2 Hessian
%
% M.Walker 11/15/2019

[nu,nv,~] = size(A);
A = reshape(A,nu*nv,3);
[Bu,dBu,ddBu] = bernsteinbasis(uvk(1),nu-1);
[Bv,dBv,ddBv] = bernsteinbasis(uvk(2),nv-1);

s = (kron(Bv,Bu)*A).';
y = x-s;
c = 0.5*sum(y.^2);

Js = [kron(Bv,dBu);kron(dBv,Bu)]*A;
grad = -Js*y;

Ay = A*y;
cu = kron(Bv,ddBu)*Ay;
cv = kron(ddBv,Bu)*Ay;
cuv = kron(dBv,dBu)*Ay;
H = Js*Js.' - [cu cuv;cuv cv];
