%% Script demonstrating surface fitting to a point cloud
%
% Plot reference surface (scaled and rotated Rosenbrock function), training and testing points, and surface
% approximation
%
% Author: Michael R. Walker II <mwalkerii@wustl.edu>
%
% This file is part of the BezBic library. See 'LICENSE' for copyright and citation info.

%% Setup
clear;

printcomment = @(varargin)fprintf('%-60s %5.1fs\n',sprintf(varargin{:}),toc);
tic;

rng(0);
rotateA = @(A,R)reshape(reshape(A,[],3)*R,size(A));

% Optionally write plot to file
output_dir = "..\doc\";
file_prefix = "plt_rosenbrock_"; 
file_prefix = "";   % no file output
plt = struct('dims',[3.4,2.8]);

% Discretization of surfaces (for visualization only)
Nuniform = 10;

%% Generate Example Surface
% Get a set of control points from a test function and generate a mesh grid of the surface

A_ref = rosenbrockTestSurface();

% Rotate control points
u = [1.09682221274760;-0.904508497187474;-0.293892626146237];
R = rotaxisangle(u/norm(u),asin(0.5*norm(u)));
A_ref = rotateA(A_ref,R);
% Get order of reference surface
nu_ref = size(A_ref,1)-1;
nv_ref = size(A_ref,2)-1;


%% Generate training and testing data
Ntr = 100;
Nte = Ntr;
sigma2 = 0.02;

uv_bound_test = 0.05;
tmp = rand(2,Ntr+Nte);
UV_tr = tmp(:,1:Ntr);
UV_te = tmp(:,Ntr+1:end)*(1-2*uv_bound_test)+uv_bound_test;


Y_tr = bezsurf.uv2x(A_ref,UV_tr) + randn(3,Ntr)*sqrt(sigma2);
wts = 0.5*(1+cos(pi*sqrt(sum((0.5-UV_tr).^2))/sqrt(1/2)));
wts = 1;
X_te = bezsurf.uv2x(A_ref,UV_te);


%% Generate reference surface
n = Nuniform;
[U,V] = meshgrid(linspace(0,1,n),linspace(0,1,n));
X_uniform = bezsurf.uv2x(A_ref,[U(:),V(:)].');

%% Fit Surface

% Regularization parameter
lambda = 1e-4;
optsA = struct('MaxIter',50,'TolX',1e-3,'TolFun',1e-5,'lambda',lambda);
optsuv = struct('MaxIter',10,'TolX',1e-2,'TolFun',1e-3,'MaxFunEvals',30,...
    'sigma',1e-1,'beta',0.5);

nu = nu_ref;
nv = nv_ref;
tic;
[Ahat,uvhat,c,exitflag,rslt] = fitsurface(Y_tr,wts,optsA,optsuv);
toc;

%% Plot fitted surface
font_sz = 8;
n = Nuniform;
[U,V] = meshgrid(linspace(0,1,n),linspace(0,1,n));
Xhat_uniform = bezsurf.uv2x(Ahat,[U(:),V(:)].');


figure('Units','inches','DefaultAxesFontSize',font_sz);
clf
dat = X_uniform.';
s1 = surf(reshape(dat(:,1),n,n),reshape(dat(:,2),n,n),reshape(dat(:,3),n,n),'EdgeColor','none','FaceColor','b','FaceAlpha',0.5);
set(gca,'FontSize',font_sz);
hold on
htr = plot3(Y_tr(1,:),Y_tr(2,:),Y_tr(3,:),'xk');
if ismac
    set(htr,'LineWidth',1);
end
hte = plot3(X_te(1,:),X_te(2,:),X_te(3,:),'.k');
dat = Xhat_uniform.';
s2 = surf(reshape(dat(:,1),n,n),reshape(dat(:,2),n,n),reshape(dat(:,3),n,n),'EdgeColor','none','FaceColor','r','FaceAlpha',0.5);
axis equal

hl = legend([htr,hte,s1,s2],'$X_\mathrm{TR}$','$S_\mathrm{TE}$','$\mathcal{S}$','$\hat{\mathcal{S}}$','Interpreter','Latex','Location','southwest');
%hl.Position(2) = hl.Position(2)+0.1;

set(gca,'View',[  215.7611   20.2290]);

plt_name = "ref_fit";
if strlength(file_prefix)
    % Squeeze the plot
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2),plt.dims]);
    % Expand axis to figure (no pad for labels)
    set(gca,'Position',[0 0.05 1 0.95])
    print('-r300',[output_dir+file_prefix+plt_name],'-depsc');
end


