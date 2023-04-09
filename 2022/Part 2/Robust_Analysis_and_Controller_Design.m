clc;
clear;
close all;
load('Assignment_Data_SC42145_2022.mat');

%% ========= Analysis of Robust Stability and Robust Performance =========
%% Construct general plant and block diagram
G = FWT(:, 1:2);        % input: beta & tau; output: omega, z.
Gd = FWT(:, 3);

% Wu and Wp
Aw = 1e-4;                          % Low-frequency disturbance attenuation
Wb = 0.3*2*pi;                      % Cut-off frequency
M = 3;                              % Infinity norm bound
s = tf('s');
Wp11 = (s/M+Wb)/(s+Wb*Aw);
Wp=[Wp11 0; 0 0.2];
Wu=[0.01 0;0 (0.005*s^2+0.0007*s+0.00005)/(s^2+0.0014*s+0.000001)];

% Wi, Wo
Wi1 = (s*(1/(16*pi))+0.2)/(s*(1/(64*pi))+1);
Wi2 = (s*(1/(16*pi))+0.2)/(s*(1/(64*pi))+1);
Wo1 = (0.05*s+0.2)/(0.01*s+1);
Wo2 = (0.05*s+0.2)/(0.01*s+1);
Wi=[Wi1 0; 0 Wi2];
Wo=[Wo1 0; 0 Wo2];

% delta_i, delta_o ---- uncertainty set
delta_i_1 = ultidyn('delta_i_1', [1 1]);
delta_i_2 = ultidyn('delta_i_2', [1 1]);
delta_o_1 = ultidyn('delta_o_1', [1 1]);
delta_o_2 = ultidyn('delta_o_2', [1 1]);
Delta_i = [delta_i_1 0; 0 delta_i_2];
Delta_o = [delta_o_1 0; 0 delta_o_2];
Delta = [Delta_i [0 0;0 0]; [0 0;0 0] Delta_o]; 

% Origianl general plant from mixed sensitivity design
P11 = [Wp; zeros(2)];
P12 = [Wp*G; Wu];
P21 = -eye(2);
P22 = -G;
P_origianl = minreal(ss([P11 P12; P21 P22]));
[K, CL, gamma, info] = hinfsyn(P_origianl, 2, 2);

%% Frequency responses of uncertainty weights
% Weights, Singular values
Wi_p = Wi * Delta_i;
Wo_p = Wo * Delta_o;

% Input weights 
figure(1)
bodemag(Wi(1, 1))
grid on;
title("Wi")

figure(2)
subplot(2, 1, 1)
bodemag(Wi_p(1, 1))     % beta
title("Beta uncertainty")
grid on;
subplot(2, 1, 2)
bodemag(Wi_p(2, 2))     % tau_e
title("Tau_e uncertainty")
grid on;

% Output weights 
figure(3)
bodemag(Wo(1, 1))
grid on;
title("Wo")

figure(4)
subplot(2, 1, 1)
bodemag(Wo_p(1, 1))     % omega
title("Omega uncertainty")
grid on;
subplot(2, 1, 2)
bodemag(Wo_p(2, 2))     % z
title("Z uncertainty")
grid on;

%Perturbed Plant
Gp = minreal((eye(2) + Wo*Delta_o) * G * (eye(2) + Wi*Delta_i));
figure(5)
sigma(Gp);
grid on;
title("Singular value of system")

%% Program the generalized plant
% General plant
%{
systemnames = 'G Wp Wu Wi Wo';
inputvar = '[ud(2); kd(2); W(2); u(2)]';
input_to_G = '[ud + u]';
input_to_Wp = '[G + kd + W]';
input_to_Wu = '[u]';
input_to_Wi = '[u]';
input_to_Wo = '[G]';
outputvar = '[Wi; Wo; Wp; Wu; -G-kd-W]';
sysoutname = 'P';
sysic;
%}

P1 = [[0 0 0 0 0 0;0 0 0 0 0 0] Wi];
P2 = [Wo*G [ 0 0 0 0;0 0 0 0] Wo*G];
P3 = [Wp*G Wp Wp Wp*G];
P4 = [zeros(2,2) [0 0 0 0;0 0 0 0] Wu];
P5 = [-G -eye(2) -eye(2) -G];
P = [P1; P2; P3; P4; P5];
P = minreal(ss(P));

% Connect uncertainty to general plant
N = minreal(lft(P, K));
M = minreal(lft(Delta, N));

%% NS, NP, RS, RP
% NS: Nominal stability
L = minreal(G * K);               % 2 outputs and 2 inputs
Det_L = minreal(1+L(1,1))*(1+L(2,2))-(L(1,2)*L(2,1));   % Det(L)
figure();
grid on;
nyquist(Det_L);

% NP: Nominal performance
omega = logspace(-3, 3, 100); % frequency space
Nf = frd(N, omega);
N_NP = Nf(5:8, 5:6);
blk = [2 4];
[mubnds, muinfo] = mussv(N_NP, blk, 'c');
muNP = mubnds(:, 1);
[muNPinf, muNPw] = norm(muNP, inf);

% RS: Robust stability
omega = logspace(-3, 3, 100); % frequency space
Nf = frd(N, omega);
N_RS = Nf(1:4, 1:4);
blk = [1 1; 1 1; 1 1; 1 1];
[mubnds, muinfo] = mussv(Nf(1:4, 1:4), blk, 'c');
muRS = mubnds(:, 1);
[muRSinf, muRSw] = norm(muRS,inf);

% RP: Robust performance
omega = logspace(-3, 3, 100); % frequency space
Nf = frd(N, omega);
blk = [1 1; 1 1; 1 1; 1 1; 2 4];
[mubnds, muinfo] = mussv(Nf, blk, 'c');
muRP = mubnds(:, 1);
[muRPinf, muRPw] = norm(muRP, inf);

figure()
loglog(muNP,muRS,muRP);
title('\mu Variations with Frequency')
grid on;
xlabel('Frequency (rad/s)');
legend('\mu_{NP}','\mu_{RS}','\mu_{RP}');
ylim([0.07, 2.5]);

%% Relation between mu for NP and H_inf norm
% Hinf obj function
Lss = minreal(series(K, G));
Ltf = tf(Lss);
S = feedback(eye(2), Ltf); 
S = minreal(S);
T = minreal(Ltf * S);       % Complementary Sensitivity (For WI)
Ti = minreal(K * S * G);    % Complementary Sensitivity (For WO)
Ob = tf(minreal([Wu*K*S; Wp*S])); 
[ninf, fpeak] = hinfnorm(Ob); %Hinfinity peak

%% ==================== Robust Controller Design ==================== 
%% DK iterations
omega = logspace(-3, 3, 100);                           % frequency space
Punc = lft(Delta, P);  
opt = dkitopt('FrequencyVector', omega,'DisplayWhileAutoIter','on');
%[K_dk, clp, dkinfo] = dksyn(Punc, 2 , 2, opt);
[K_dk, clp, dkinfo] = musyn(Punc, 2 , 2);
K_dk = minreal(ss(K_dk));

%% Check NS, NP, RS, and RP using SSV
% Nominal Stability
N_dk = minreal(lft(P, K_dk));
eig_dk = max(real(eig(N_dk)));

% Nominal performance
Nf_dk=frd(N_dk,omega);
blk = [2 4];
[mubnds, muinfo] = mussv(Nf_dk(5:8, 5:6),blk,'c');
muNP_dk = mubnds(:, 1);
[muNPinf_dk, muNPw_dk] = norm(muNP_dk,inf);

% Robust Stability
blk = [1 1; 1 1; 1 1; 1 1];
[mubnds, muinfo] = mussv(Nf_dk(1:4, 1:4), blk, 'c');
muRS_dk = mubnds(:, 1);
[muRSinf_dk,muRSw_dk] = norm(muRS_dk,inf);

% Robust Performance
blk = [1 1; 1 1; 1 1; 1 1; 2 4];
[mubnds, muinfo] = mussv(Nf_dk(1:8, 1:6), blk, 'c');
muRP_dk = mubnds(:, 1);
[muRPinf_dk, muRPw_dk] = norm(muRP_dk, inf);

% Normalized SSV bounds for NP, RS, RP
figure()
loglog(muNP_dk, muRS_dk, muRP_dk);
title('\mu')
xlabel('Frequency (rad/s)');
legend('\mu_{NP}-DK', '\mu_{RS}-DK', '\mu_{RP}-DK');
ylim([0.07, 2.5]);

%% Time-domain simulation scenarios (MIMO vs. DK)
% DK controller
L = G * K_dk;                                       % Open loop system 
S = minreal(feedback(eye(2), L));                   % Sensitivity function
CL = minreal(feedback(L, eye(2)));                  % Close loop tf

% Mixed-Sensitivity controller
Ln = G * K;                                         % Open loop system
S2 = minreal(feedback(eye(2), Ln));                 % Sensitivity function
CL2 = feedback(Ln, eye(2));                         % Close loop tf

% ==== Time domain plots
% Disturbance Rejection
figure()
step(S2 * Gd);
hold on;
step(S * Gd);
hold off;
grid on;
title('Disturbance Rejection')
legend('MS', 'DK')
xlabel('Time');
xlim([0, 200]);

% Reference Tracking
figure()
step(CL2(:,1));
hold on;
step(CL(:,1));
hold off;
grid on;
title('Reference Tracking')
legend('MS', 'DK')
xlabel('Time');
xlim([0, 200]);

% ==== Frequency domain plots
figure()
MS_plant = G * K;
DK_plant = G * K_dk;
bodemag(MS_plant(:, 1));
hold on;
bodemag(DK_plant(:, 1));
hold off;
grid on;
legend('MS', 'DK')
title('Bode diagram of open loop system');

loopmixed = loopsens(G, K);
loopdk = loopsens(G, K_dk);
% Senstivity function 
S_ms = minreal(inv(eye(2) + MS_plant));
S_dk = minreal(inv(eye(2) + DK_plant));
% Complementary sensitivity
T_ms = minreal(MS_plant * (inv(eye(2) + MS_plant)));
T_dk = minreal(DK_plant * (inv(eye(2) + DK_plant)));

% Visualizaion
figure();
bodemag(S_ms(:, 1));
hold on;
bodemag(S_dk(:, 1));
hold off;
grid on;
legend('MS', 'DK');
title('Sensitivity');

figure();
bodemag(T_ms(:, 1));
hold on;
bodemag(T_dk(:, 1));
hold off;
grid on;
legend('MS', 'DK');
title('Complementary Sensitivity');