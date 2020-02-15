function A = rosenbrockTestSurface(a,b,sf)
% ROSENBROCKTESTSURFACE - Generate Bezier surface control points approximating the Rosenbrock function
%
% The Rosenbrock function is defined
%  $f(x,y) = (a-x)^2 + b\left(y-x^2\right)^2$
%
% We sample f(x,y) over a mesh grid of points,  (x,y) \in S, to obtain tuples (x,y,f(x,y)), or (x,y,z). We select (u,v)
% coordinates for each point by scaling and translating (x,y). A least-squares solution fitting the control points, A,
% is available in closed form.
%
% A = rosenbrockTestSurface()
% A = rosenbrockTestSurface(a,b)
% A = rosenbrockTestSurface(a,b,sf)
%  Optionally specify a scaling to the z-coordinate (DEFAULT = 1/625)
%
% M.Walker 12/3/2019

% Polynomial order
ordr = 4;

% Number of points to generate for fitting. This should be large enough to ensure B is full-rank.
n = 100;

if nargin < 3
    sf = 1/625;
    if nargin < 2
        b = 100;
        if nargin < 1
            a = 1;
        end
    end
end

% Rosenbrock function
% Also, scale the coordinates to accomodate rotations
[X,Y] = meshgrid(linspace(-2,2,n),linspace(-1,3,n));
Z = (b*(Y-X.^2).^2 + (a-X).^2)*sf;

% Generate u,v mesh by scaling x,y pairs
U = X - min(X(:));
U = U./max(U(:));
V = Y - min(Y(:));
V = V./max(V(:));

% Fit controlpoints to scaled x,y pairs
Bu = bernsteinbasis(U(:),ordr);
Bv = bernsteinbasis(V(:),ordr);
% kron(Bv,Bu) ROWS only
B = reshape(Bu.*permute(Bv,[1 3 2]),n^2,(ordr+1)^2);
A = B\[X(:),Y(:),Z(:)];

A = reshape(A,ordr+1,ordr+1,3);