classdef UKF < handle
    % UKF Implementation for vehicle state estimation
    properties
    end

    methods (Access = public)
    end

    methods (Access = private)
        
    end

    methods (Static)
        function [Xi,W] = SigmaPoints(xm,P,kappa)
            % Compute SigmaPoints and respective weights
            n = numel(xm);
            Xi = zeros(n,2*n+1);
            W = zeros(2*n+1,1);
            Xi(:,1) = xm;
            W(1) = kappa/(n+kappa);
            
            U = chol((n+kappa)*P);

            for k=1:n
                Xi(:,k+1) = xm + U(k,:)';
                W(k+1) = 1/(2*(n+kappa));
            end

            for k=1:n
                Xi(:,n+k+1) = xm - U(k,:)';
                W(n+k+1) = 1/(2*n+kappa);
            end
        end

        function [xm,xcov] = UT(Xi,W)
            [n,kmax] = size(Xi);
            xm = 0;
            for k=1:kmax
                xm = xm + W(k) * Xi(:,k);
            end

            xcov = zeros(n,n);
            for k=1:kmax
                xcov = xcov + W(k) * (Xi(:,k) - xm) * (Xi(:,k) - xm)';
            end
            
        end
    end
end