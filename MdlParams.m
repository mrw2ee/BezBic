classdef MdlParams
    % MDLPARAMS - Model parameters
    %   Retain all model parameters
    %
    % MdlParams(A,uv)
    %  A  : (nu+1)x(nv+1)x3
    %  uv : 2 x nx
    % MdlParams(A,uv,sigma2)
    %  sigma2 : scalar
    %
    % MdlParams(z,nu,nv)
    %  z = [vect(A);vect(uv)] or [vect(A);vect(uv);sigma2]
    %
    % Default sigma2 : NaN
    %
    % M.Walker 12/2/2019
    
    properties
        A;
        uv;
        sigma2 = NaN;
    end
    properties (Dependent)
        nu;
        nv;
        nx;
        d;
        z;
    end
    
    methods
        function obj = MdlParams(A,uv,sigma2)
            if isscalar(uv)
                if ~isscalar(sigma2)
                    error('Input dimensions inconsistent');
                end
                uv0 = (uv+1)*(sigma2+1)*3+1;
                obj.A = reshape(A(1:uv0-1),uv+1,sigma2+1,3);
                if mod(numel(A)-numel(obj.A),2) == 0
                    % Sigma not included
                    obj.uv = reshape(A(uv0:end),2,[]);
                    obj.sigma2 = NaN;
                else
                    obj.uv = reshape(A(uv0:end-1),2,[]);
                    obj.sigma2 = A(end);
                end
            else
                obj.A = A;
                obj.uv = uv;
                if nargin > 2
                    obj.sigma2 = sigma2;
                end
            end
        end
            
        function x = get.nu(obj)
            x = size(obj.A,1)-1;
        end
        function x = get.nv(obj)
            x = size(obj.A,2)-1;
        end
        function x = get.nx(obj)
            x = size(obj.uv,2);
        end
        function x = get.d(obj)
            x = numel(obj.A) + numel(obj.uv) + 1;
        end
        function x = get.z(obj)
            if isnan(obj.sigma2)
                x = [obj.A(:);obj.uv(:)];
            else
                x = [obj.A(:);obj.uv(:);obj.sigma2];
            end
        end
    end
end

