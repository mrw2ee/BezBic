%% bernsteinbasisTest
%
classdef bezsurfTest < matlab.unittest.TestCase
    properties (TestParameter)
        order = struct('small', 2, 'medium', 4,'large', 6);
    end
    methods (Test)
        function uvcostGrad(testCase,order)
            % Check first derivative of Bernstein polynomials
            L = 1e2;
            n = order;
            %n = 4;
            u = linspace(0,1,L).';
            [Y,Yd] = bernsteinbasis(u,n);
            
            % Compute derivative numerically
            Neps = 50;
            eps_all = 10.^linspace(-10,-1,Neps);
            err_all = zeros(L,n+1,numel(eps_all));
            for itt = 1:Neps
                eps = eps_all(itt);
                d = sign(randn(L,1));
                Yprim = bernsteinbasis(u+eps*d,n);
                err_all(:,:,itt) = (Yprim-Y)/eps - Yd.*d;
            end
            err_all = permute(err_all,[3,2,1]);
            
            figure
            loglog(eps_all,abs(reshape(err_all,Neps,[])))
            ylabel('Error')
            xlabel('Step size: \epsilon')
            title(sprintf('Gradient of Bernstein basis vectors: Order %d',n))
            grid on
            
            % Check linear growth for reasonable eps
            idx = eps_all < 1e-4 & eps_all > 1e-6;
            rslt = log(abs(reshape(err_all(idx,:,:),sum(idx),[]))./eps_all(idx).');
            % At some values u the gradient is nearly exact. Log error will not grow linearly with log eps. Throw out these values
            idx = max(rslt) < 0 & min(rslt) < -3;
            rslt = var(rslt(:,~idx));
            fprintf('1st Derivative: Order %d, Max error %e, Ignoring %d of %d\n',order,max(rslt),sum(idx),numel(idx));
            verifyLessThan(testCase,rslt,1e-3);
            
        end
        
        function testDerivative2(testCase,order)
            % Check second derivative of Bernstein polynomials
            L = 1e2;
            n = order;
            u = linspace(0,1,L).';
            [~,Yd,Ydd] = bernsteinbasis(u,n);
            
            % Compute derivative numerically
            Neps = 50;
            eps_all = 10.^linspace(-10,-1,Neps);
            err_all = zeros(L,n+1,numel(eps_all));
            for itt = 1:Neps
                eps = eps_all(itt);
                d = sign(randn(L,1));
                [~,Yprim] = bernsteinbasis(u+eps*d,n);
                err_all(:,:,itt) = (Yprim-Yd)/eps - Ydd.*d;
            end
            err_all = permute(err_all,[3,2,1]);
            
            figure
            loglog(eps_all,abs(reshape(err_all,Neps,[])))
            ylabel('Error')
            xlabel('Step size: \epsilon')
            title(sprintf('2nd derivative of Bernstein basis vectors: Order %d',n))
            grid on
            
            % Check linear growth for reasonable eps
            idx = eps_all < 1e-4 & eps_all > 1e-6;
            rslt = log(abs(reshape(err_all(idx,:,:),sum(idx),[]))./eps_all(idx).');
            % At some values u the gradient is nearly exact. Log error will not grow linearly with log eps. Throw out these values
            idx = max(rslt) < 0 & min(rslt) < -3;
            rslt = var(rslt(:,~idx));
            fprintf('2nd Derivative: Order %d, Max error %e, Ignoring %d of %d\n',order,max(rslt),sum(idx),numel(idx));
            if sum(~idx)>0
                verifyLessThan(testCase,rslt,1e-3);
            end
        end
        
        
        
        function testBezierSurface(testCase)
            % Simple sanity check
            L = 20;
            n = 3;
            
            [X,Y] = meshgrid(0:n,0:n);
            Z = sin(X/2*pi).*sin(Y/2*pi);
            
            A = X;
            A(:,:,2) = Y;
            A(:,:,3) = Z;
            
            Y = bernsteinbasis(linspace(0,1,L).',n);
            rslt = zeros(L,L,3);
            for it = 1:3
                rslt(:,:,it) = Y*A(:,:,it)*Y.';
            end
            rslt = permute(rslt,[2,1,3]);
            Z = reshape(rslt,[],3).';
            
            [U,V] = meshgrid(linspace(0,1,L).',linspace(0,1,L).');
            Zhat = bezsurf.uv2x(A,[U(:),V(:)].');
            
            verifyEqual(testCase,Zhat,Z,'AbsTol',2e-5);
        end
    end
end