%% Simple Vehicle Model based FGO Simulation
% 
% Step 0: EKF, UKF based state estimation + RTS Smoother
% Step 1: Test for Offline Batch Optimization
% Step 2: Test for Online Optimization (Sliding Window)
% 
% 2D Vehicle Model: [x, y, v, psi]
% 
% [State Propagation Eqn]
% x_(k+1) = x_k + v_k * cos(psi_k) * dt
% y_(k+1) = y_k + v_k * sin(psi_k) * dt
% v_(k+1) = v_k + a_k * dt
% psi_(k+1) = psi_k + psidot * dt

%% Read data
clear; close all; clc;
load('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\VehicleLocalizationandMapping\output.mat');
load('D:\2021_Spring\research\3Secondz\Dataset\OnlineLDR\VehicleLocalizationandMapping\lane.mat');

% Scenario 1 DSTART = 14800; DEND = 16000 --> 기준은 lane에 있는 t
% Scenario 2 DSTART = 17000; DEND = 21800

% Alter timestamp value for output.ubloxgps.t(79)
predT = output.ubloxgps.t(78) + 0.05;
delta = predT - output.ubloxgps.t(79);
output.ubloxgps.t(79:end) = output.ubloxgps.t(79:end) + delta;
% plot(output.ubloxgps.t(1:100) - output.ubloxgps.t(1));

%% Test EKF with RTS Smoother
sol = struct();
sol.EKF = EKF(output,1,230000);
sol.EKF.optimize();
sol.EKF.plotRes();

%% Test UKF with RTS Smoother

%% Test Offline Batch Optimization

%% Test Online Sliding Window Optimization