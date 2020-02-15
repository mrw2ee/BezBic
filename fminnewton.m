function [xnext,varargout] = fminnewton(fun,x0,opts)
% fminnewton - Newton minimization with Armijo backtracking
%
% x = fminnewton(fun,x0)
% x = fminnewton(fun,x0,opts)
%   Opts - structure with 'MaxIter','TolX','TolFun','MaxFunEvals','sigma','beta'
%   'MaxIter' - Maximum nubmer of Newton iterations
%   'MaxFunEvals' - Maximum number of function evals (may exceed iterations due to Armijo backtracking)
%   'TolX','TolFun' - Tolerance on parameters and cost
%   'sigma' - Threshold scaling [1e-5,1e-1].
%   'beta' - Stepsize reduction [0.1,0.5]
% [x, fx, exitflag, output] = fminnewton(fun,x0,opts)
%   
%
% M.Walker 11/15/2019
if nargin < 3
    opts = struct('MaxIter',100*numel(x0),'TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',100*numel(x0),...
        'sigma',1e-1,'beta',0.5);
end
sigma = opts.sigma;
beta = opts.beta;

exitflag = 0;
xnext = x0;
[fnext,g,H] = fun(xnext);
funevals = 1;
for iteration = 1:opts.MaxIter
    % Preserve start at iteration
    xlast = xnext;
    flast = fnext;
    
    % Newton direction:
    d = -H\g;
    graddotd = d.'*g;
    if d.'*g > 0
        % Not descent direction, switch to steepest descent
        graddotd = 0.5;
        d = -graddotd*g;
    end
    
    xnext = xlast + d;
    [fnext,g,H] = fun(xnext);
    funevals = funevals + 1;
    
    % Armijo backtracking
    bm = 1;
    while (flast-fnext)<-sigma*bm*graddotd
        bm = bm*beta;
        xnext = xlast + bm*d;
        [fnext,g,H] = fun(xnext); 
        funevals = funevals + 1;
        
        if norm(xnext-xlast,'fro') < opts.TolX ||  funevals > opts.MaxFunEvals
            % Give up!, but ensure fnext is best
            if fnext > flast
                xnext = xlast;
                fnext = flast;
            end
            break;
        end
    end
    
    % Are we done?
    if (norm(xnext-xlast,'fro') < opts.TolX || abs(fnext-flast) < opts.TolFun)
        % Converged!
        exitflag = 1;
        break;
    elseif funevals > opts.MaxFunEvals
        % Give up!
        break;
    end
end

varargout = cell(1,nargout-1);
if nargout > 1
    if nargout > 2
        if nargout > 3
            varargout{3} = struct('iterations',iteration,'funcCount',funevals);
        end
        varargout{2} = exitflag;
    end
    varargout{1} = fnext;
end