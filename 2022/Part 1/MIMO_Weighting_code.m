clc;
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
% Check given wind data
figure();
plot(WindData(:, 1), WindData(:, 2));
grid on;
title("Wind data in time domain");
xlabel("time (s)");
ylabel("magnitude");

figure();
% show two peaks
fs = 100;               % sampling frequency
n = length(WindData);
X = fft(WindData(:, 2));
f = (0:n-1)*(fs/n);     %frequency range
power = abs(X).^2/n;    %power
plot(f,power)
xlim([-5,120])

% show lower frequency and higher frequency
figure();
Fs = 100;                % Sampling frequency             
Y = fft(WindData(:, 2));
L = length(WindData);
f = (0:length(Y)-1)*Fs/length(Y);
P1 = abs(Y);
subplot(2,2,1);
plot(f,P1) ;
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("f (Hz)");
ylabel("|P1(f)|");
xlim([0 10]);
title("Wind data in frequency domain");
subplot(2,2,2);
plot(f,P1) ;
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("f (Hz)");
ylabel("|P1(f)|");
xlim([0 2]);
title("Wind data in frequency domain");
subplot(2,2,3);
plot(f,P1) ;
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("f (Hz)");
ylabel("|P1(f)|");
xlim([0 0.5]);
title("Wind data in frequency domain");
subplot(2,2,4);
plot(f,P1) ;
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("f (Hz)");
ylabel("|P1(f)|");
xlim([0 0.3]);
title("Wind data in frequency domain");

% For output discurbance rejection W_p, supporting theory can be
% seen from book page 94-95
s = tf([1 0], 1);
Mi = 2;
Ai = 1e-3;
omega_B = 0.08;                 % Chose this so that S is bounded
Wp = (s/Mi + omega_B) / (s + omega_B * Ai);
figure();
bodemag(Wp);                    % Frequency respond of Mp
yline(0,'--');
grid on;
title("Bode diagram of Wp");

% For controller W_u
% one higher order filter and one lower order filter
Wu_beta = makeweight(0.25, 0.5*6.28, 100);       % cut off frequency 0.5Hz accordiong to text
Wu_tau_e = makeweight(50, 0.5*6.28, 0.25); 
Wu = [tf(Wu_beta) 0; 0 tf(Wu_tau_e)]; 
figure()
bodemag(Wu_beta, Wu_tau_e);
yline(0,'--');
grid on;
title("Bode diagram of Wu");
legend("Beta", "Tau_e", "0 dB");

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
bodemag(S);
hold on;
grid on;
bodemag(T);
title("Sensitivity and Complementary Sensitivity");
legend("Sensitivity","Complementary sensitivity");
hold off;

Sd = Gd / (1 + L);
KSd = tf(K) * Gd / (1 + L);
KS = K * S;

% Sensitivity to Wp
%   S is bounded by 1/wp
figure();
bodemag(S);
%bodemag(Sd);
hold on;
grid on;
bodemag(1/Wp);
hold off;
legend("Sensitivity", "1/Wp");
title("Sensitivity to weight --- Wp");

% Sensitivity to Wu
%   KS is bounded by 1/wu
figure();
bodemag(KS(1));
%bodemag(KSd(1));
hold on;
grid on;
bodemag(1/Wu(1, 1));
hold off;
legend("Sensitivity", "1/Wu (beta)");
title("Sensitivity to weight --- Wu(beta)");

figure();
bodemag(KS(1));
%bodemag(KSd(2));
hold on;
grid on;
bodemag(1/Wu(2, 2));
hold off;
legend("Sensitivity", "1/Wu (tau)");
title("Sensitivity to weight --- Wu(tau)");

% Controller to Wu beta (aiming low frequency input)
figure()
bodemag(K(1))           % beta control
hold on;
grid on;
bodemag(Wu(1, 1));
hold off;
legend("Controller", "Wu (beta)");
title("Controller to weight --- Wu(beta)");

% Controller to Wu tau (aiming high frequency input)
figure()
bodemag(K(2))           % tau control
hold on;
grid on;
bodemag(Wu(2, 2));
hold off;
legend("Controller", "Wu (tau)");
title("Controller to weight --- Wu(tau)");

%% System simulation
figure()
n = length(WindData);
t = linspace(0,10*60,n);
u = WindData(:, 2);
lsim(inv(1+L)*Gd, u, t);
yline(0,'--');
title("Wind data input");

figure()
% Sin signal for testing
texting_signal = sin((t*2*pi)/1000);        % Ampitude2, frequency 2e-3pi
lsim(inv(1+L)*Gd, texting_signal, t);
yline(0,'--');
title("Sine signal input");