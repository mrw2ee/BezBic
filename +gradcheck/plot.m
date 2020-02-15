function plot(x0,fctn,varargin)
% plot - Plot gradient error vs step size
%
% plot(x,fctn)
%
% x - Initial starting point (may be multi-dimensional vector)
% fctn - Function handle returning the function value and gradient
%   [fx,dfx] = fctn(x)
%
% Parameter-value pairs:
%
% 'eps' - Specify vector of step sizes
% 'N' - Specify number of random directions to evaluate (DEFAULT = 5)
% 'nonnegative' - TF, Constrain x>= 0. DEFAULT = false
%   Note, this will preclude consistent direction vectors and eps values. Results are corrected internally.
% 'hessian' - T/F, Evaluate Hessian. DEFAULT = false
%
% See also gradcheck.singledirection
%
% M.Walker 11/17/2019

default_eps = 10.^linspace(-10,2,50);
default_N = 5;

p = inputParser;
addRequired(p,'x',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
addRequired(p,'fctn',@(x)validateattributes(x,{'function_handle'},{'nonempty'}));
addParameter(p,'eps',default_eps,@(x)validateattributes(x,{'numeric'},{'row','positive'}))
addParameter(p,'N',default_N,@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'}));
addParameter(p,'nonnegative',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'hessian',false,@(x)validateattributes(x,{'logical'},{'scalar'}));

parse(p,x0,fctn,varargin{:});

Ncolors = size(get(gca,'ColorOrder'),1);

eps = p.Results.eps;
Nlines = 2*(1+p.Results.hessian);
for didx = 1:p.Results.N
    [err,eps_actual] = gradcheck.singledirection(x0,fctn,eps,[],p.Results.nonnegative,p.Results.hessian);
    % Plot
    loglog(eps_actual,err);
    hold on
    set(gca,'ColorOrderIndex',mod(Ncolors-Nlines-1+get(gca,'ColorOrderIndex'),Ncolors)+1);
    if didx == 1
        xlabel('\epsilon')
        ylabel('Error');
        grid on;
    end
    drawnow;
end
set(gca,'ColorOrderIndex',mod(Ncolors+Nlines-1+get(gca,'ColorOrderIndex'),Ncolors)+1);