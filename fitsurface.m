function [A,uv,varargout] = fitsurface(X,wts,optsA,optsuv)
% fitsurface - Fit Bezier surface and determine model order automatically
%
% [A,uv] = fitsurface(X,wts)
%   X - 3xN data
%   wts - weights Nx1 (or scalar)
%   A - Matrix of control points
%   uv - 2xN matrix of surface coordinates for the data
%
% Optional inputs
% fitsurface(..., optsA, optsuv)
%   optsA - structure of parameters for A update: MaxIter, TolX, TolFun, lambda (Tikhonov regularization)
%   optsuv - structure of parameters for uv update. See fminnewton
%
% Optional additional outputs:
% [A,uv,sigma2,exitflag,log]
%   sigma2 - estimated noise variance
%   exitflag - Positive result indicates converged
%   log - Structure with vectored history of sigma2, uv update set size, and size of A.
%
% M.Walker 12/17/2019

if nargin < 3 || isempty(optsA)
    optsA = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-2,'lambda',1e-2);
end
if nargin < 4 || isempty(optsuv)
    optsuv = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-3,'MaxFunEvals',30,...
        'sigma',1e-1,'beta',0.5);
end

nu = 1;
nv = 1;

lambda = optsA.lambda;
Xmu = median(X,2);
X = X - Xmu;
uv = bezsurf.x2uvsvd(X.*wts);
[A,sigma2] = ModelSelect(X,wts,uv,nv,nu,lambda);

exitflag = 0;
sigma2_hist = zeros(optsA.MaxIter,1);
step_hist = sigma2_hist;
Asize = zeros(optsA.MaxIter+1,1);
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
    [A,sigma2] = ModelSelect(X,wts,uv,size(A,2)-1,size(A,1)-1,lambda);
    
    Asize(iteration) = size(A,1)*size(A,2);
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
            Asize(iteration+1:end) = [];
            varargout{3} = struct('sigma2',sigma2_hist,'step',step_hist,'Asize',Asize);
        end
        varargout{2} = exitflag;
    end
    varargout{1} = sigma2;
end
end
function [A,sigma2] = updateA(X,wts,Bv,Bu,lambda)
% UPDATEA
%
% X : 3xN
% wts : 1xN
% Bv : (n_v +1) x N
% Bu : (n_u +1) x N
% lambda : scalar
nx = size(Bu,2);
B = reshape(permute(Bu,[1 3 2]).*permute(Bv,[3 1 2]),[],nx);
Bwtssqrd = B.*(wts.^2);
A = (Bwtssqrd*B.' + lambda*eye(size(B,1)))\(Bwtssqrd*X.');
sigma2 = norm((A.'*B - X).*wts,'fro')^2/(3*size(X,2));
end

function [A,sigma2] = ModelSelect(X,wts,uv,nv,nu,lambda)
% ModelSelect
Bu = bernsteinbasis(uv(1,:).',nu).';
Bup = bernsteinbasis(uv(1,:).',nu+1).';
Bv = bernsteinbasis(uv(2,:).',nv).';
Bvp = bernsteinbasis(uv(2,:).',nv+1).';

[A,sigma2] = updateA(X,wts,Bv,Bu,lambda);
[Au,sigma2u] = updateA(X,wts,Bv,Bup,lambda);
[Av,sigma2v] = updateA(X,wts,Bvp,Bu,lambda);
[Auv,sigma2uv] = updateA(X,wts,Bvp,Bup,lambda);

nx = size(X,2);
d = 1 + 2*nx + 3*(nu+1+[0,1,0,1]).*(nv+1+[0,0,1,1]);
bic =- 3*nx*log([sigma2,sigma2u,sigma2v,sigma2uv]) - d*log(nx);
[~,idx] = max(bic);
switch(idx)
    case 2
        A = Au;
        nu = nu+1;
        sigma2 = sigma2u;
    case 3
        A = Av;
        nv = nv+1;
        sigma2 = sigma2v;
    case 4
        A = Auv;
        nu = nu+1;
        nv = nv+1;
        sigma2 = sigma2uv;
end
A = reshape(A,nu+1,nv+1,3);
end
