classdef EKF < handle
    %EKF Implementation for vehicle state estimation

    properties
        t
        accel
        vel
        gyro
        yaw
        gps
        Xkk
        Xkkp1 = []
        Pkk
        Pkkp1 = []
        Fk = []
        Hk = []
        idx1
        idx2
        Q = diag([0.05^2, 0.05^2, 0.03^2, 1e-4])
        R = diag([1e-4, 1e-4])
    end

    methods (Access = public)
        %% Constructor
        function obj = EKF(output,idx1,idx2)            
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
            % OPTIMIZE Run EKF throughout data points
            obj.Xkk = zeros(4,obj.idx2-obj.idx1+1);
            obj.Xkk(:,1) = [obj.gps(1);obj.gps(2);obj.vel(1);obj.yaw(1)];
            Xkk_ = obj.Xkk(:,1);
            Pkk_ = diag([1e-4,1e-4,1e-4,1e-6]);
            obj.Pkk = zeros(16,obj.idx2-obj.idx1+1);
            obj.Pkk(:,1) = reshape(Pkk_,[],1);
            
            obj.Xkkp1 = zeros(4,obj.idx2-obj.idx1);
            obj.Pkkp1 = zeros(16,obj.idx2-obj.idx1);
            obj.Fk = zeros(16,obj.idx2-obj.idx1);
            obj.Hk = zeros(8,obj.idx2-obj.idx1);
            for i=1:obj.idx2-obj.idx1
                [Xkkp1_,Pkkp1_,Fk_] = obj.predict(Xkk_,Pkk_,obj.Q,obj.accel(i),obj.gyro(i),obj.t(i+1) - obj.t(i));
                [Xkk_,Pkk_,Hk_] = obj.update(Xkkp1_,Pkkp1_,obj.R,obj.gps(:,i));
                obj.Xkkp1(:,i) = Xkkp1_;
                obj.Pkkp1(:,i) = reshape(Pkkp1_,[],1);
                obj.Xkk(:,i+1) = Xkk_;
                obj.Pkk(:,i+1) = reshape(Pkk_,[],1);
                obj.Fk(:,i) = reshape(Fk_,[],1);
                obj.Hk(:,i) = reshape(Hk_,[],1);
            end
        end
        
        %% Plot Results
        function plotRes(obj)
            %PLOTRES Plot state estimation results
            
            % X-Y
            figure(1); hold on; grid on; axis equal;
            plot(obj.Xkk(1,:),obj.Xkk(2,:),'k.'); 
            plot(obj.gps(1,:),obj.gps(2,:),'r.');
            xlabel('X'); ylabel('Y');

            % Yaw Comparison
            figure(2); hold on; grid on;
            plot(obj.t - obj.t(1),180/pi * obj.Xkk(4,:),'k.');
            plot(obj.t - obj.t(1),180/pi * obj.yaw,'r.');
            xlabel('Time(s)'); ylabel('Yaw(deg)');            
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

        %% Update    
        function [Xkk,Pkk,Hk] = update(Xkkp1,Pkkp1,Rk,Zk)
            Hk = [1, 0, 0, 0;
                  0, 1, 0, 0];
            yk = Zk - Hk * Xkkp1;
            Sk = Hk * Pkkp1 * Hk' + Rk;
            Kk = Pkkp1 * Hk' / Sk;
            Xkk = Xkkp1 + Kk * yk;
            Pkk = (eye(4) - Kk * Hk) * Pkkp1;
        end

    end
end