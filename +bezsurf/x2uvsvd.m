function UV = x2uvsvd(X,tol)
% x2uvsvd - Map points to uv coordinates using SVD
%
% UV = x2uvsvd(X,tol)
% 
% X - 3xN input
% tol - Optional offset (default = 0)
%
% UV - 2xN output
%
% Where UV \in [tol,1-tol].
%
% M.Walker 4/23/2019

if nargin < 2
    tol = 0;
end

% Center data prior to SVD
Xmu = median(X,2);
X = X - Xmu;
[~,~,UV] = svds(X,2);
UV = UV.';
UV = UV - min(UV,[],2);
UV = UV.*((1-2*tol)./max(UV,[],2)) + tol;