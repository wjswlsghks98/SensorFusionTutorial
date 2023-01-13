classdef UKF < handle
    % UKF Implementation for vehicle state estimation
    properties
        t
        accel
        vel
        gyro
        yaw
        gps
        idx1
        idx2
        Xkk
        Pkk
        kappa = -1
        Q = diag([0.05^2, 0.05^2, 0.03^2, 1e-4])
        R = diag([1e-4, 1e-4])
    end

    methods (Access = public)
        function obj = UKF(output,idx1,idx2)
            obj.t = output.rtimu.t(idx1:idx2);
            obj.accel = [output.rtimu.ax(idx1:idx2); output.rtimu.ay(idx1:idx2)];
            obj.vel = [output.rtvel.vx(idx1:idx2); output.rtvel.vy(idx1:idx2)];
            obj.gyro = output.rtimu.wz(idx1:idx2);
            obj.yaw = output.rtimu.rz(idx1:idx2);
            obj.gps = [output.rtgps.x(idx1:idx2); output.rtgps.y(idx1:idx2)];
            obj.idx1 = idx1; obj.idx2 = idx2;
            % Initialize
            obj.initialize();
        end
        
        %% Optimize
        function obj = optimize(obj)
            obj.Xkk = zeros(4,obj.idx2-obj.idx1+1);
            obj.Xkk(:,1) = [obj.gps(1);obj.gps(2);obj.vel(1);obj.yaw(1)];
            Xkk_ = obj.Xkk(:,1);
            Pkk_ = diag([1e-4,1e-4,1e-4,1e-6]);
            obj.Pkk = zeros(16,obj.idx2-obj.idx1+1);
            obj.Pkk(:,1) = reshape(Pkk_,[],1);
            
            for i=1:obj.idx2-obj.idx1
                [Xi,W] = obj.SigmaPoints(Xkk_,Pkk_,obj.kappa);
                
            end
        end
    end

    methods (Access = private)
        %% Pre-Process input data
        function obj = initialize(obj)
            new_accel = obj.accel(1,:) .* cos(obj.yaw) + obj.accel(2,:) .* sin(obj.yaw);
            obj.accel = new_accel;
            new_vel = obj.vel(1,:) .* cos(obj.yaw) + obj.vel(2,:) .* sin(obj.yaw);
            obj.vel = new_vel;
        end
    end

    methods (Static)
        %% Predict
        function [Xkkp1,Pkkp1,Fk] = predict(Xkk,Pkk,Qk,ak,wk,dt)
            vk = Xkk(3); psik = Xkk(4);
            Xkkp1 = Xkk + [vk * cos(psik) * dt;
                           vk * sin(psik) * dt;
                           ak * dt;
                           wk * dt];
            Fk = [1, 0, cos(psik), -vk * sin(psik);
                  0, 1, sin(psik), vk * cos(psik);
                  0, 0, 1, 0;
                  0, 0, 0, 1];
            Pkkp1 = Fk * Pkk * Fk' + Qk;            
        end
        
        %%
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
        
        %%
        function [xm,xcov] = UT(Xi,W)
            % Unscented Transform
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