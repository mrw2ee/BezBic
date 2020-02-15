%% testuvsolver
%
function tests = testuvsolver
    tests = functiontests(localfunctions);
end

function testUvOptimization(testCase)

% Generate points on known surface
n = 20;
A = rosenbrockTestSurface();
pad = 0.1;
[U,V] = meshgrid(linspace(pad,1-pad,n),linspace(pad,1-pad,n));
X = bezsurf.uv2x(A,[U(:),V(:)].');

% Add noise to surface coordinates
rng(0);
sigma2 = 0.1^2;
UV0 = [U(:),V(:)] + randn(numel(U),2)*sqrt(sigma2);
UVstar = [U(:),V(:)].';
% Check Gradient & Hessian
opts = struct('MaxIter',100,'TolX',1e-6,'TolFun',1e-8,'MaxFunEvals',100,...
        'sigma',1e-1,'beta',0.5);
for k = 1:n
    uv = UV0(k,:).';
    x = X(:,k);
    uv_hat = fminnewton(@(uv)point2surfacel2(uv,x,A),uv,opts);
    
    verifyEqual(testCase,uv_hat,UVstar(:,k),'AbsTol',1e-5);
end
end