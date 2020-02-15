%% Script printing table of performance metrics
%
%
% Author: Michael R. Walker II <mwalkerii@wustl.edu>
%
% This file is part of the BezBic library. See 'LICENSE' for copyright and citation info.

%% Setup
clear;

pctic = tic;
printcomment = @(varargin)fprintf('%-60s %5.1fs\n',sprintf(varargin{:}),toc(pctic));

rng(0);

% Regularization parameter
lambda = 1e-7;
% Number of training and testing points
Ntr = 100;
Nte = Ntr;

uv_bound_test = 0.05;

Nit = 100;
Ntr_all = [50,100,500,1000];
eps_all = [1e-2,5e-2,1e-1,5e-1].^2;

if exist('Paper_Metrics.mat','file')
    load('Paper_Metrics','T');
else
    T = table('Size',[0,8],'VariableTypes',{'double','double','double','double','double','double','double','double'},...
        'VariableNames',{'Ntr','sigma2_Y','Nit','Asize','sigma2hat_tr','sigma2hat_te','sigma2_tr','t'});
    for nidx = 1:numel(Ntr_all)
        printcomment('Nidx %d of %d',nidx,numel(Ntr_all));
        Ntr = Ntr_all(nidx);
        for eidx = numel(eps_all):-1:1
            eps = eps_all(eidx);
            for it = 1:Nit
                %%
                optsA = struct('MaxIter',10,'TolX',1e-5,'TolFun',1e-5,'lambda',lambda);
                optsuv = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',20,...
                    'sigma',1e-1,'beta',0.1);
                
                [X_tr,S_te,A_ref] = generatepoints(Ntr,Nte,eps,uv_bound_test);
                % Training noise: roughly eps/3
                S_tr = footpoints(A_ref,X_tr);
                sigma2_tr = norm(S_tr - X_tr,'fro')^2/(3*Ntr);
                
                fittic = tic;
                [Ahat,~,c,~,rslt] = fitsurface(X_tr,1,optsA,optsuv);
                t = toc(fittic);
                
                % Testing error: on order of c
                Shat_te = footpoints(Ahat,S_te);
                sigma2hat_te = norm(Shat_te - S_te,'fro')^2/(3*Nte);
                
                T = [T;{Ntr,eps,numel(rslt.sigma2),rslt.Asize(end),c,sigma2hat_te,sigma2_tr,t}];
            end
        end
    end
    save('Paper_Metrics','T');
end

%%
Tmb = table2means(T);

%%
dat = Tmb{:,[1:6,end]};
dat(:,end) = dat(:,end)*1e3;
input = struct;
input.data = dat;
input.tableColumnAlignment = '|cc|ccccc|';
input.dataFormat = {'%g',1,'%.1E',1,'%g',2,'%.1E',2,'%.1f',1};
input.tableColLabels = {'$n_\textrm{TR}$', '$\sigma^2_Y$','iter.','size','$\hat{\sigma}^2_\textrm{TR}$','$\hat{\sigma}^2_\textrm{TE}$','ms'};
input.tableBorders = 0;
input.tableLabel = 'itrslts';
input.tableCaption = 'Effects of Problem Size and Noise on Results';

fprintf('\n\n');
if ~exist('latexTable','file')
    disp(Tmb);
else
    tmp = latexTable(input);
    
    npreamble = 4;
    nfooter = 3;
    ndata = numel(tmp)-npreamble-nfooter;
    hlinesat = npreamble + [0,1,1+(4:4:ndata+1)];
    for itt = numel(hlinesat):-1:1
        tmp = [tmp(1:hlinesat(itt));{'\hline'};tmp(hlinesat(itt)+1:end)];
    end
    npreamble = npreamble+1;
    % Add spacing around header
    tmp{npreamble} = [tmp{npreamble},'\rule{0pt}{2.4ex}'];
    tmp{npreamble+1} = [tmp{npreamble+1},'[2pt]'];
    
    disp(char(tmp))
end
%%
function [Y_tr,Y_te,A_ref] = generatepoints(Ntr,Nte,sigma2,uv_bound_test)
rotateA = @(A,R)reshape(reshape(A,[],3)*R,size(A));
A_rosenbrock = rosenbrockTestSurface();
u = randn(3,1);
u = u/norm(u);
x = rand(1);
R = rotaxisangle(u,x,'uniform');
A_ref = rotateA(A_rosenbrock,R);

tmp = rand(2,Ntr+Nte);
UV_tr = tmp(:,1:Ntr);
UV_te = tmp(:,Ntr+1:end)*(1-2*uv_bound_test)+uv_bound_test;

Y_tr = bezsurf.uv2x(A_ref,UV_tr) + randn(3,Ntr)*sqrt(sigma2);
Y_te = bezsurf.uv2x(A_ref,UV_te);
end

function [X,varargout] = meshsurf(A,n)
% MESHSURF - Mesh Bezier surface with control points A
if nargin < 2
    n = 10;
end
[U,V] = meshgrid(linspace(0,1,n),linspace(0,1,n));
X = bezsurf.uv2x(A,[U(:),V(:)].');
if nargout>1
    varargout = {[U(:),V(:)].'};
end
end

function [s1,s2,htr,hte] = plotsurfs(Ahat,A_ref,Y_tr,Y_te)
n = 10;
Xhat = meshsurf(Ahat,n);
Xref = meshsurf(A_ref,n);

dat = Xref.';
s1 = surf(reshape(dat(:,1),n,n),reshape(dat(:,2),n,n),reshape(dat(:,3),n,n),'EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
hold on
dat = Xhat.';
s2 = surf(reshape(dat(:,1),n,n),reshape(dat(:,2),n,n),reshape(dat(:,3),n,n),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);
if nargin > 2
    htr = plot3(Y_tr(1,:),Y_tr(2,:),Y_tr(3,:),'xk');
    hte = plot3(Y_te(1,:),Y_te(2,:),Y_te(3,:),'.k');
end
end

function Yhat = footpoints(Ahat,Y_te)

optsuv = struct('MaxIter',50,'TolX',1e-3,'TolFun',1e-3,'MaxFunEvals',60,...
    'sigma',1e-1,'beta',0.5);
Nte = size(Y_te,2);
n = 10;
[X_rough,UV_rough] = meshsurf(Ahat,n);
xidx = zeros(1,Nte);
parfor n = 1:Nte
    [~,xidx(n)] = min(sum((X_rough - Y_te(:,n)).^2));
end
UV0 = UV_rough(:,xidx);

% Update UV
parfor n = 1:Nte
    xi = Y_te(:,n);
    uvi = UV0(:,n);
    uv(:,n) = fminnewton(@(uv)point2surfacel2(uv,xi,Ahat),uvi,optsuv);
end

% Performance metric
Yhat = bezsurf.uv2x(Ahat,uv);
end

function Tout = table2means(Tin)

Tout = Tin;
Tout(1:end,:) = [];
Ntr_all = unique(Tin.Ntr);
eps_all = unique(Tin.sigma2_Y);
for nidx = 1:numel(Ntr_all)
    Ntr = Ntr_all(nidx);
    for eidx = numel(eps_all):-1:1
        eps = eps_all(eidx);
        idx = Tin{:,1} == Ntr & Tin{:,2} == eps;
        %dat = [dat;mean(T{idx,[1:6,8]})];
        Tout = [Tout;num2cell(mean(Tin{idx,:}))];
        
    end
end
end