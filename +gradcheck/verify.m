function verify(testCase, x0,fctn,varargin)
% verify - Unit test case for gradient (and Hessian)
%
% verify(testCase, x, fctn)
%
% testCase - If empty, creates TestCase.forInteractiveUse
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
% 'TolSlope' - Abs. tolerance error on slope (DEFAULT = 1e-1)
% 'Tol0' - Abs. tolerance on initial error (DEFAULT = 1e-3). Note, multiply by 20 for Hessian.
%
% See also gradcheck.plot, gradcheck.singledirection
%
% M.Walker 11/17/2019

default_eps = 10.^linspace(-10,2,50);
default_N = 5;
default_TolSlope = 1e-1;
default_Tol0 = 1e-3;

p = inputParser;
addRequired(p,'x',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
addRequired(p,'fctn',@(x)validateattributes(x,{'function_handle'},{'nonempty'}));
addParameter(p,'eps',default_eps,@(x)validateattributes(x,{'numeric'},{'row','positive'}))
addParameter(p,'N',default_N,@(x)validateattributes(x,{'numeric'},{'scalar','positive','integer'}));
addParameter(p,'nonnegative',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'hessian',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'TolSlope',default_TolSlope,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'Tol0',default_Tol0,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));

parse(p,x0,fctn,varargin{:});
if isempty(testCase)
    testCase = matlab.unittest.TestCase.forInteractiveUse;
end

Nlines = 2*(1+p.Results.hessian);
eps = p.Results.eps;
err = zeros(numel(eps),Nlines,p.Results.N);
eps_actual = zeros(numel(eps),p.Results.N);
for didx = 1:p.Results.N
    [err(:,:,didx),eps_actual(:,didx)] = gradcheck.singledirection(...
        x0,fctn,eps,[],p.Results.nonnegative,p.Results.hessian);
end

% Quantify errors
dat = diff(log10(err))./diff(log10(permute(eps_actual,[1,3,2])));

% Average across N samples since single outliers are likely

% Check first row is small.
testCase.verifyLessThan(abs(mean(dat,3)-1),p.Results.TolSlope);
testCase.verifyLessThan(mean(err(1,1:0.5*end,:),3),p.Results.Tol0);
testCase.verifyLessThan(mean(err(1,0.5*end+1:end,:),3),p.Results.Tol0*20);

