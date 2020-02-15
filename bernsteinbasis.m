function [B,varargout] = bernsteinbasis(u,n)
% bernsteinbasis - Evaluate the Bernstein basis at requested points
%
% Syntax
%
% B = bernsteinbasis(u,n)
%
% For p-element column vector u, returns px(n+1) matrix of sampled basis vectors. Depending on u these vectors are
% typically independent but not orthogonal.
%
% [B,Bd] = bernsteinbasis(u,n)
% [B,Bd,Bdd] = bernsteinbasis(u,n)
%
% Optionally return first and second derivatives of the basis w.r.t. u
%
% M.Walker 4/10/2019

validateattributes(u,{'numeric'},{'column'},1)
if nargin < 2
    n = 2;
end
validateattributes(n,{'numeric'},{'positive','integer','scalar'},2)

ui = u.^(0:n);
ui2 = (1-u).^(n:-1:0);
coef = factorial(n)./(factorial(0:n).*factorial(n:-1:0));
B = coef.*(ui.*ui2);

L = numel(u);

% Note, naive computation of the derivative may include inverse factors of u. However, these will all be multiplied by a
% 0 coefficient. To avoid NaNs, we simply ommit these terms.
if nargout > 1
    if nargout > 2
        % B"
        varargout{2} = [zeros(L,2),coef(3:end).*(2:n).*(1:n-1).*ui(:,1:end-2).*ui2(:,3:end)] ...
            +[zeros(L,1),-2*coef(2:end-1).*(1:n-1).*(n-1:-1:1).*ui(:,1:end-2).*ui2(:,3:end),zeros(L,1)]...
            + [coef(1:end-2).*(n:-1:2).*(n-1:-1:1).*ui(:,1:end-2).*ui2(:,3:end),zeros(L,2)];
    end
    % B'
    varargout{1} = [zeros(L,1),coef(2:end).*(1:n).*ui(:,1:end-1).*ui2(:,2:end)] ...
        - [coef(1:end-1).*(n:-1:1).*ui(:,1:end-1).*ui2(:,2:end),zeros(L,1)];
end