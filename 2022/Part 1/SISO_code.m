clc;
clear;
close all;
load('Assignment_Data_SC42145_2022.mat')

% ============ SISO Analysis and Control Design ============
% construct model
%   since we want tocontrol beta to increase w_r, system
%   is modeified accordingly
A1 = A;
B1 = B(:, 1);       % beta - w_r
C1 = -1* C(1, :);   % *-1 according to hint
D1 = 0;
[b, a] = ss2tf(A1, B1, C1, D1);
ol_sys = tf(b, a);
sys_zpk = zpk(ol_sys);

%% Q1 Open loop analysis
zeros_ol = zero(sys_zpk);
poles_ol = pole(sys_zpk);
[gm, pm ,wcg, wcp] = margin(ol_sys);
%plot_system_all(ol_sys);

%% Q2 Project requirements
cl_system = feedback(ol_sys, 1);
%plot_system_all(cl_system);     

%% Q3 Design controller for closed-loop system
Kp = 1;           % 1.1
Ki = 0;           % 0.5
Kd = 0;
PID_controller = pid(Kp,Ki,0,0);
controlled_ol_sys = ol_sys * PID_controller;

% Find the threshold of Kp 
for kp = [5, 10, 15, 16, 20]
    buffer = pid(kp,0,0,0) * ol_sys;
    stepinfo(feedback(buffer, 1));
    % system become unstable when kp reaches 16
end

% Find the threshold of Ki 
clear buffer;
for ki = [0.2,0.25,0.3,0.5,0.75,1,2]
    buffer = pid(1,ki,0,0) * ol_sys;
    %figure();
    %step(feedback(buffer, 1));
    stepinfo(feedback(buffer, 1));
    % system become unstable when ki reaches 1
end

%% Design the combine PI controller
clc;
clear buffer;
controlled_ol_system = pid(1, 0.25, 0, 0) * ol_sys;
controlled_cl_system = feedback(controlled_ol_system, 1);
plot_system_all(controlled_cl_system)

%% Q4 Rejecting output disturbance
clear buffer;
A2 = A;
B2 = B(:, 3);       % disturbance v -> w_r
C2 = C(1, :);   
D2 = 0;
[b2, a2] = ss2tf(A2, B2, C2, D2);
ol_sys2 = tf(b2, a2);
sys_zpk2 = zpk(ol_sys2);

% new system don't need any -1 since the formation changed
A3 = A;
B3 = B(:, 1);       % beta - w_r
C3 = C(1, :);       % No need for -1 anymore
D3 = 0;
[b3, a3] = ss2tf(A3, B3, C3, D3);
ol_sys3 = tf(b3, a3);

% Basic information of the system
zeros_ol2 = zero(sys_zpk2);
poles_ol2 = pole(sys_zpk2);

% Construct system with distrubuence

figure();
for p = [0.5, 0.75, 1, 1.25, 1.5, 2, 5]
    combined_system_dis = ol_sys2 / (1 - pid(p, 0.25, 0, 0) * ol_sys3);
    step(combined_system_dis);
    hold on;
    grid on;
end
title("Step respond of different p");
legend("p=0.5", "p=0.75", "p=1", "p=1..25", "p=1.5", "p=2", "p=5");


figure();
for k = [0.25, 0.5, 0.75, 1, 2]
    combined_system_dis = ol_sys2 / (1 - pid(2, k, 0, 0) * ol_sys3);
    step(combined_system_dis);
    hold on;
    grid on;
end
title("Step respond of different k");
legend("k=0.25", "k=0.5", "k=0.75", "k=1", "k=2");

% plot the optimal system
optimal_combined_system_dis = ol_sys2 / (1 - pid(2, 1, 0, 0) * ol_sys3);
plot_system_all(optimal_combined_system_dis)

%% Visualization function definition
function [Gain_Margin, Phase_Margin, Wcg, Wcp] = plot_system_all(sys)
figure()            % create new figure
% Step respond
subplot(2,2,1);
time = 0:0.01:500;
stepplot(sys, time);                               
grid on;
stepinfo(sys)

% Bode diagram with margin
subplot(2,2,2);
margin(sys);
grid on;
[Gain_Margin,Phase_Margin, Wcg, Wcp] = margin(sys);

% Nyquist plot 
subplot(2,2,3)
nyquist(sys);
grid on;

% Root locus 
subplot(2,2,4)
rlocus(sys);
grid on; 
end