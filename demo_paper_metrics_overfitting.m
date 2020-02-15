%% Scipt generating plots contrasting order selection with fixed-order models
%
% Author: Michael R. Walker II <mwalkerii@wustl.edu>
%
% This file is part of the BezBic library. See 'LICENSE' for copyright and citation info.


%% Setup
clear;
pctic = tic;
printcomment = @(varargin)fprintf('%-60s %5.1fs\n',sprintf(varargin{:}),toc(pctic));

rng(0);

output_dir = "..\doc\";
file_prefix = "plt_fitting_";
file_prefix = "";
plt = struct('dims',[1.6,2.8]);
font_sz = 8;
if exist('Paper_MetricsOverUnder.mat','file')
    load('Paper_MetricsOverUnder','Tbic','Tfixed');
else

    % Regularization parameter
    lambda = 1e-7;
    % Number of training and testing points
    Ntr = 100;
    Nte = Ntr;
    % Sampling noise
    eps = 0.5^2;
    uv_bound_test = 0.05;

    Ntr_all = [100];
    eps_all = 10.^(-3:0.5:0);
    Nitt = 10;

    Tbic = table('Size',[0,9],'VariableTypes',{'double','double','string','double','double','double','double','double','double'},...
        'VariableNames',{'Ntr','sigma2_Y','surf','Nit','Asize','sigma2hat_tr','sigma2hat_te','sigma2_tr','t'});
    Tfixed = arrayfun(@(x)Tbic,1:5,'UniformOutput',false);

    optsA = struct('MaxIter',10,'TolX',1e-5,'TolFun',1e-5,'lambda',lambda);
                optsuv = struct('MaxIter',30,'TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',10,...
                    'sigma',1e-1,'beta',0.1);


    for nidx = 1:numel(Ntr_all)
        Ntr = Ntr_all(nidx);
        for eidx = numel(eps_all):-1:1
            eps = eps_all(eidx);
            printcomment('Eidx %d of %d',eidx,numel(eps_all));
            for it = 1:Nitt
                %%
                for this_surf = {'plane','rosen'}
                    if isequal(this_surf{1},'plane')
                        [X_tr,S_te,A_ref] = generatepointsPlane(Ntr,Nte,eps,uv_bound_test);
                    else
                        [X_tr,S_te,A_ref] = generatepointsRosen(Ntr,Nte,eps,uv_bound_test);
                    end
                    % Training noise: roughly eps/3
                    S_tr = footpoints(A_ref,X_tr);
                    sigma2_tr = norm(S_tr - X_tr,'fro')^2/(3*Ntr);

                    fittic = tic;
                    [Ahat,~,c,~,rslt] = fitsurface(X_tr,1,optsA,optsuv);
                    t = toc(fittic);

                    % Testing error: on order of c
                    Shat_te = footpoints(Ahat,S_te);
                    sigma2hat_te = norm(Shat_te - S_te,'fro')^2/(3*Nte);

                    Tbic = [Tbic;{Ntr,eps,this_surf{1},numel(rslt.sigma2),rslt.Asize(end),c,sigma2hat_te,sigma2_tr,t}];

                    for odr = 1:numel(Tfixed)
                        fittic = tic;
                        [Ahat,~,c,~,rslt] = fitsurfacefixedorder(X_tr,1,odr,odr,optsA,optsuv);
                        t = toc(fittic);
                        % Testing error: on order of c
                        Shat_te = footpoints(Ahat,S_te);
                        sigma2hat_te = norm(Shat_te - S_te,'fro')^2/(3*Nte);
                        Tfixed{odr} = [Tfixed{odr};{Ntr,eps,this_surf{1},numel(rslt.sigma2),numel(Ahat(:,:,1)),c,sigma2hat_te,sigma2_tr,t}];
                    end
                end
            end
        end
    end
    save('Paper_MetricsOverUnder','Tbic','Tfixed');
end

%%
Tmb = table2means(Tbic);
Tmf = cellfun(@table2means,Tfixed,'UniformOutput',false);

%%
fhdlp = figure(1);
clf;
s = "plane";
Ntr = 100;
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
Tall = [Tmf,{Tmb}];
for tidx = 1:numel(Tall)
    T = Tall{tidx};
    idx = T.surf == s & T.Ntr == Ntr;
    loglog(T.sigma2_Y(idx),T.sigma2hat_te(idx),[markers{mod(tidx-1,numel(markers))+1},'-']);
    hold on;
end

grid on
xlabel('$\sigma^2_Y$','Interpreter','Latex');
ylabel('$\hat{\sigma}^2_\textrm{TE}$','Interpreter','Latex');
ahdlp = gca;



fhdlr = figure(2);
clf;
s = "rosen";

markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
Tall = [Tmf,{Tmb}];
for tidx = 1:numel(Tall)
    T = Tall{tidx};
    idx = T.surf == s & T.Ntr == Ntr;
    loglog(T.sigma2_Y(idx),T.sigma2hat_te(idx),[markers{mod(tidx-1,numel(markers))+1},'-']);
    hold on;
end


xlabel('$\sigma^2_Y$','Interpreter','Latex');
ylabel('$\hat{\sigma}^2_\textrm{TE}$','Interpreter','Latex');
legend([arrayfun(@(x)sprintf('%d c-pts',Tfixed{x}.Asize(1)),1:numel(Tfixed),'UniformOutput',false),'ours'],'Location','SouthEast');
%title(sprintf('%s, %d points',s,Ntr))

if strlength(file_prefix)
    %ylim(yl);
    ahdlr = gca;

    for fhdl = [fhdlp,fhdlr]
        fhdl.Units = 'inches';
        pos = fhdl.Position;
        set(fhdl,'Position',[pos(1:2),plt.dims]);
    end
    linkaxes([ahdlp,ahdlr],'xy');
    for ahdl = [ahdlp,ahdlr]
        set(ahdl,'FontSize',font_sz);
    end
    grid on

    plt_name = "plane";
    print(fhdlp,'-r300',[output_dir+file_prefix+plt_name],'-depsc');

    plt_name = "rosen";
    print(fhdlr,'-r300',[output_dir+file_prefix+plt_name],'-depsc');
end

%%
function [Y_tr,Y_te,A_ref] = generatepointsRosen(Ntr,Nte,sigma2,uv_bound_test)
rotateA = @(A,R)reshape(reshape(A,[],3)*R,size(A));
A = rosenbrockTestSurface();
u = randn(3,1);
u = u/norm(u);
x = rand(1);
R = rotaxisangle(u,x,'uniform');
A_ref = rotateA(A,R);

tmp = rand(2,Ntr+Nte);
UV_tr = tmp(:,1:Ntr);
UV_te = tmp(:,Ntr+1:end)*(1-2*uv_bound_test)+uv_bound_test;

Y_tr = bezsurf.uv2x(A_ref,UV_tr) + randn(3,Ntr)*sqrt(sigma2);
Y_te = bezsurf.uv2x(A_ref,UV_te);
end

function [Y_tr,Y_te,A_ref] = generatepointsPlane(Ntr,Nte,sigma2,uv_bound_test)
rotateA = @(A,R)reshape(reshape(A,[],3)*R,size(A));
A = zeros(2,2,3);
A(:,:,1) = [-2 -2;2 2];
A(:,:,2) = [-1 3;-1 3];

u = randn(3,1);
u = u/norm(u);
x = rand(1);
R = rotaxisangle(u,x,'uniform');
A_ref = rotateA(A,R);

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

Tout = table('Size',[0,size(Tin,2)+1],'VariableTypes',[arrayfun(@(x)class(Tin{1,x}),1:size(Tin,2),'UniformOutput',false),'double'],...
    'VariableNames',[Tin.Properties.VariableNames,'sigma2hat_te_stdv']);


Ntr_all = unique(Tin.Ntr);
eps_all = unique(Tin.sigma2_Y);
surf_all = unique(Tin.surf);
for nidx = 1:numel(Ntr_all)
    Ntr = Ntr_all(nidx);
    for eidx = numel(eps_all):-1:1
        eps = eps_all(eidx);
        for sidx = 1:numel(surf_all)
            surf = surf_all(sidx);
            idx = Tin{:,1} == Ntr & Tin{:,2} == eps & Tin{:,3} == surf;
            %dat = [dat;mean(T{idx,[1:6,8]})];
            %Tout = [Tout;num2cell(mean(Tin{idx,1:2})),{surf},num2cell(mean(Tin{idx,4:end}))];
            Tout = [Tout;num2cell(mean(Tin{idx,1:2})),{surf},num2cell(mean(Tin{idx,4:end})),sqrt(var(Tin{idx,end-2}))];
        end
        
    end
end
end