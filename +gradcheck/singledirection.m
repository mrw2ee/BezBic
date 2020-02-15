function [err,eps_actual] = singledirection(x0,fctn,eps,d,nonnegative,evalhess)
%% graderr - Quantify graident accuracy
%
% SYNTAX 
%
% [err,eps] = graderr(x,fctn,d,eps)
%
% x - Initial starting point (may be multi-dimensional vector)
% fctn - Function handle returning the function value and gradient
%   [fx,dfx] = fctn(x)
%
% eps - Specify step sizes
%
%  [err,eps] = graderr(x,fctn,d,eps,NONNEGATIVE) - Constrain x>= 0. DEFAULT = false
%   Note, this will preclude consistent direction vectors and eps values. Results are corrected internally. Returned eps
%   represents actual step size.
%
% graderr(...,NONNEGATIVE,EVALHESSIAN) - Additionally evaluate Hessian
%   If true, we assume fctn supports three outputs
%   [fx,dfx,H] = fctn(x)
%
%
% ALGORITHM
%
% To quantify error we use the Taylor expansion
%
% $$ f( x + d \epsilon) \approx f(x) + \nabla f(x) ^T d $$
% 
% This motivates the error function
%
% $$ e(\epsilon) := \left| \frac{f( x + d \epsilon) - f(x)}{\epsilon} - \nabla f(x)^T d\right| $$
% 
% Since this requires evaluation of the function at $f( x + d \epsilon)$ at all $\epsilon$, we also compute the gradient
% $\nabla f( x + d \epsilon)$ which we can use to form an additional error function. Due to the symmetry of the aboslute
% value, we have
%
% $$ e(\epsilon) := \left| \frac{f( x + d \epsilon) - f(x)}{\epsilon} - \nabla f(x + d \epsilon)^T d\right| $$
%
% Which we report as a second column of the output error.
%
% To quantify the Hessian we consider the Taylor expansion of the Gradient
%
% $$ \nabla f( x + d \epsilon) \approx \nabla f(x) + \epsilon \nabla^2 f(x) d $$
%
% This motivates the error function
%
% $$ e(\epsilon) := \left\| \frac{\nabla f( x + d \epsilon) - \nabla f(x)}{\epsilon} - \nabla^2 f(x) d\right\|_{\ell 2} $$
%
% Here we replace the absolute values with the vector norm. Again, due to symmetry, we append two columns to the
% returned error matrix.
%
% See also gradcheck.plot
%
% M.Walker 11/20/2019

% Evaluate initial cost and gradient
if ~evalhess
    % Gradient only
    err = zeros(numel(eps),2);
    [f0,df0] = fctn(x0);
else
    % Evaluate Hessian
    err = zeros(numel(eps),4);
    [f0,df0,H0] = fctn(x0);
end

eps_actual = eps;

if isempty(d)
    d = randn(size(x0));
    d = d/norm(d,'fro');
end
% For each stepsize...
for eidx = 1:numel(eps)
    % Evaluate cost and gradient
    x = x0 + eps(eidx)*d;
    if nonnegative && (min(x,[],'all')<0)
        x = max(x,0);
        % Note, direction and eps changed...
        this_d = x - x0;
        this_eps = norm(this_d,'fro');
        this_d = this_d/this_eps;
    else
        this_d = d;
        this_eps = eps(eidx);
    end
    if ~evalhess
        [f,df] = fctn(x);
    else
        [f,df,H] = fctn(x);
    end

    % Quantify accuracy
    eps_actual(eidx) = this_eps;
    % Gradient
    tmp = (f-f0)/this_eps;
    err(eidx,1) = abs(tmp - sum(df0.*this_d,'all'));
    err(eidx,2) = abs(tmp - sum(df.*this_d,'all'));
    if evalhess
        % Hessian
        tmp = (df-df0)./this_eps;
        err(eidx,3) = norm(tmp - H0*this_d,'fro');
        err(eidx,4) = norm(tmp - H*this_d,'fro');
    end
end