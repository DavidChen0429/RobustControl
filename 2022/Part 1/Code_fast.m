clear;
close all;
load('Assignment_Data_SC42145_2022.mat');

% =============== MIMO Weighting Design ===============
%% Construct the generalized plant
A_MIMO = A;
B_MIMO = B(:, 1:3);               % Input: Beta, Tau_e, V
C_MIMO = C(1, :);                 % Output: omega_r
D_MIMO = [0 0 0];

% Construct state space and transfer functions
[num1, den1] = ss2tf(A_MIMO, B_MIMO, C_MIMO, D_MIMO, 1); % Control input: Beta
[num2, den2] = ss2tf(A_MIMO, B_MIMO, C_MIMO, D_MIMO, 2); % Control input: Tau_e
[num3, den3] = ss2tf(A_MIMO, B_MIMO, C_MIMO, D_MIMO, 3); % Exogenous input: V

V_wr_tf = tf(num3, den3);               
Beta_wr_tf = tf(num1, den1);                        % Control input 1
Tau_e_wr_tf = tf(num2, den2);                       % Control input 2
Gd = V_wr_tf;                                       % Gd
Input_to_output_tf = [Beta_wr_tf Tau_e_wr_tf];      % G

%% Design Weights
% For output discurbance rejection W_p, supporting theory can be
% seen from book page 94-95
s = tf([1 0], 1);
Mi = 2;
Ai = 1e-3;
omega_B = 0.08;                 % Chose this so that S is bounded
Wp = (s/Mi + omega_B) / (s + omega_B * Ai);
% For controller W_u
% one higher order filter and one lower order filter
Wu_beta = makeweight(0.25, 0.5*6.28, 100);       % cut off frequency 0.5Hz accordiong to text
Wu_tau_e = makeweight(100, 0.5*6.28, 0.25); 
Wu = [tf(Wu_beta) 0; 0 tf(Wu_tau_e)]; 

%% Synthesize controllers with weights
P = minreal(augw(Input_to_output_tf, Wp, Wu, [])); % Construct generalized plant
[K, CL, gamma] = hinfsyn(P, 1, 2);

% Construct open system
G = Input_to_output_tf;
%Gd has already assigned
K = minreal(K);
L = G*K;

% Senitivity function 
% This plot shows that S and T both have expected quality for low and high
% frequency
S = 1 / (1 + L);
T = L / (1 + L);
Sd = Gd / (1 + L);
KSd = tf(K) * Gd / (1 + L);
KS = K * S;

%% Simulate the system
figure()
n = length(WindData);
t = linspace(0,10*60,n);
u = WindData(:, 2);
lsim(inv(1+L)*Gd, u, t);
title("Wind data input");

figure()
% Sin signal for testing
texting_signal = sin((t*2*pi)/1000);        % Ampitude2, frequency 2e-3pi
lsim(inv(1+L)*Gd, texting_signal, t);
title("Sine signal input");