function [A,uv,varargout] = fitsurfacefixedorder(X,wts,nu,nv,optsA,optsuv)
% fitsurfacefixedorder - Fit Bezier surface with fixed model order
%
% [A,uv] = fitsurfacefixedorder(X,wts,nu,nv)
%   X - 3xN data
%   wts - weights Nx1 (or scalar)
%   nu,nv - Surface order in two dimensions (positive integers)
%   A - Matrix of control points
%   uv - 2xN matrix of surface coordinates for the data
%
% Optional inputs
% fitsurface(..., optsA, optsuv)
%   optsA - structure of parameters for A update: MaxIter, TolX, TolFun, lambda (Tikhonov regularization)
%   optsuv - structure of parameters for uv update. See fminnewton
%
% M.Walker 12/17/2019
if nargin < 5 || isempty(optsA)
    optsA = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-2,'lambda',1e-2);
end
if nargin < 6 || isempty(optsuv)
    optsuv = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-3,'MaxFunEvals',30,...
        'sigma',1e-1,'beta',0.5);
end

lambda = optsA.lambda;
Xmu = median(X,2);
X = X - Xmu;
uv = bezsurf.x2uvsvd(X.*wts);
[A,sigma2] = updateA(X,uv,nu,nv,wts,lambda);

exitflag = 0;
sigma2_hist = zeros(optsA.MaxIter,1);
step_hist = sigma2_hist;
for iteration = 1:optsA.MaxIter
    cost_last = sigma2;
    uv_last = uv;
    A_last = A;
    
    % Update UV
    parfor n = 1:size(X,2)
        xi = X(:,n);
        uvi = uv_last(:,n);
        uv(:,n) = fminnewton(@(uv)point2surfacel2(uv,xi,A),uvi,optsuv);
    end
    
    % Update A
    [A,sigma2] = updateA(X,uv,nu,nv,wts,lambda);
    
    sigma2_hist(iteration) = sigma2;
    step_hist(iteration) = norm(uv-uv_last,'fro');
    
    % Should we stop?
    if sigma2 > cost_last
        % Cost increased!
        uv = uv_last;
        A = A_last;
        sigma2 = cost_last;
        break;
    elseif (cost_last-sigma2 < optsA.TolFun) 
        % Converged!
        exitflag = 1;
        break;
    elseif step_hist(iteration)<optsA.TolX
        % Converged!
        exitflag = 2;
        break;
    end
end
A = A + permute(Xmu,[2 3 1]);

varargout = cell(1,nargout-1);
if nargout > 2
    if nargout > 3
        if nargout > 4
            sigma2_hist(iteration+1:end) = [];
            step_hist(iteration+1:end) = [];
            varargout{3} = struct('sigma2',sigma2_hist,'step',step_hist);
        end
        varargout{2} = exitflag;
    end
    varargout{1} = sigma2;
end
end
function [A,sigma2] = updateA(X,uv,nu,nv,wts,lambda)
% UPDATEAt
%
% X : 3xN
% uv : 2 x nx matrix
% nu,nv : basis order (positive integers)
% wts : 1xN
% lambda : scalar
B = computeB(uv,nu,nv);
Bwtssqrd = B.*(wts.^2);
A = (Bwtssqrd*B.' + lambda*eye(size(B,1)))\(Bwtssqrd*X.');
sigma2 = norm((A.'*B - X).*wts,'fro')^2/(3*size(X,2));
A = reshape(A,nu+1,nv+1,3);
end
