clc;
clear all;
close all;
load('Assignment_Data_SC42145.mat');
%% 3.1

s = tf('s');
A1=A;
B1=B(:,1:2);
C1=C;
D1=0;
ltisys=ss(A1,B1,C1,D1);
G=ss(ltisys);

% Wp, Wu
wB=0.4*2*pi;
DR=1/10000;
M=1.8;
Wp=[(s/M +wB)/(s+wB*DR) 0; 0 0.2];
Wu=[0.01 0;0 (0.005*s^2+0.0007*s+0.00005)/(s^2+0.0014*s+0.000001)];

%Generalized plant for mixed-sensitivity design
P11 = [Wp; zeros(2)];
P12 = [Wp*G; Wu];
P21 = -eye(2);
P22 = -G;
P = minreal(ss([P11 P12; P21 P22]),0.0001);
[K,CL,gamma,info] = hinfsyn(P,2,2);
K=minreal(balreal(K)); %mixed-sensitivity controller

% Wi, Wo
Wi=[(s*(1/(16*pi))+0.3)/(s*(1/(64*pi))+1) 0; 0 (s*(1/(16*pi))+0.3)/(s*(1/(64*pi))+1)];
Wo=[(0.05*s+0.2)/(0.01*s+1) 0;0 (0.05*s+0.2)/(0.01*s+1)];

% delta_i, delta_o
delta_i_1 = ultidyn('delta_i_1',[1 1]);
delta_i_2 = ultidyn('delta_i_2',[1 1]);
delta_o_1 = ultidyn('delta_o_1',[1 1]);
delta_o_2 = ultidyn('delta_o_2',[1 1]);
Delta_i=[delta_i_1 0;0 delta_i_2];
Delta_o=[delta_o_1 0;0 delta_o_2];
Delta=[Delta_i [0 0;0 0];[0 0;0 0] Delta_o]; %Perturbation block Delta

% Definition of the new generalized plant
P1=[[0 0 0 0 0 0;0 0 0 0 0 0] Wi];
P2=[Wo*G [ 0 0 0 0;0 0 0 0] Wo*G];
P3=[Wp*G Wp Wp Wp*G];
P4=[zeros(2,2) [0 0 0 0;0 0 0 0] Wu];
P5=[-G -eye(2) -eye(2) -G];
P=[P1;P2;P3;P4;P5];
P=minreal(ss(P),0.0001);

N=lft(P,K); %N
N=minreal(N,0.00001);
M=lft(Delta, N); %M
M=minreal(M,0.00001);


% Weights, Singular values
Wi_p=Wi*Delta_i;
Wo_p=Wo*Delta_o;
figure(1)
bodemag(Wi(1,1))
figure(2)
subplot(2, 1, 1)
bodemag(Wi_p(1,1))
subplot(2,1,2)
bodemag(Wi_p(2,2))

figure(3)
bodemag(Wo(1,1))
figure(4)
subplot(2, 1, 1)
bodemag(Wo_p(1,1))
subplot(2,1,2)
bodemag(Wo_p(2,2))

%Perturbed Plant
Gp = (eye(2)+Wo*Delta_o)*G*(eye(2)+Wi*Delta_i);
figure(5)
sigma(Gp);

%Nominal Stability
OL=minreal(G*K,0.00001);
n=minreal((1+OL(1,1))*(1+OL(2,2))-(OL(1,2)*OL(2,1)),0.00001);
pol=pole(OL); %Checking number of OL RHP poles P
zer=tzero(OL); %Checking Number of CL RHP poles Z=Number of OL RHP zeros
figure(6)
nyquist(n)
NS = max(real(eig(N)));

% Nominal performance
omega=logspace(-3,3,100); %frequency range
Nf=frd(N,omega); %Convert to frequency domain
blk=[2 4]; %Size of perturbation block
[mubnds,muinfo]=mussv(Nf(5:8,5:6),blk,'c'); %get u bounds
muNP=mubnds(:,1);
[muNPinf,muNPw]=norm(muNP,inf); %normalize and find peak
figure(7)
loglog(muNP);
grid on;
title('\mu NP')
xlabel('Frequency (rad/s)');
ylim([0.07,2.5]);

% Robust Stability
blk=[1 1;1 1;1 1;1 1];
[mubnds,muinfo]=mussv(Nf(1:4,1:4),blk,'c');
muRS=mubnds(:,1);
[muRSinf,muRSw]=norm(muRS,inf);
figure(8)
loglog(muRS);
grid on;
title('\mu RS')
xlabel('Frequency (rad/s)');
ylim([0.07,2.5]);

%Robust Performance
blk=[1 1; 1 1; 1 1; 1 1; 2 4];
[mubnds,muinfo]=mussv(Nf(1:8,1:6),blk,'c');
muRP=mubnds(:,1);
[muRPinf,muRPw]=norm(muRP,inf);
figure(9)
loglog(muRP);
grid on;
title('\mu RP')
xlabel('Frequency (rad/s)');
ylim([0.07,2.5]);

figure(10)
loglog(muNP,muRS,muRP);
title('\mu Variations with Frequency')
grid on;
xlabel('Frequency (rad/s)');
legend('\mu_{NP}','\mu_{RS}','\mu_{RP}');
ylim([0.07,2.5]);

% Hinf obj function
Lss = minreal(series(K,G),0.00001);
Ltf = tf(Lss); %Open loop Tf
S = feedback(eye(2),Ltf); %Sensitvity
S = minreal(S);
T=minreal(Ltf*S); %Complementary Sensitivity(For WI)
Ti=minreal(K*S*G); %Complementary Sensitivity(For WO)
Ob = [Wu*K*S; Wp*S]; 
Ob = minreal(Ob, 0.00001);
Ob = tf(Ob);
[ninf,fpeak] = hinfnorm(Ob); %Hinfinity peak

bound=frd(ninf*ones(1, 100),omega);
figure(11)
semilogx(muNP, bound)
legend('\mu_{NP}', 'H_{\infty} peak');
ylim([0,2]);
grid on;
title('\mu Synthesis Nominal Performance vs H_{\infty} norm peak')
xlabel('Frequency (rad/s)');

%% 3.2

%D-K Iterations
P_dk=lft(Delta,P);  
P_dk=minreal(P_dk);
figure(12)
opt=dkitopt('FrequencyVector', omega,'DisplayWhileAutoIter','on');
[K_dk,CL_dk,dkinfo] = dksyn(P_dk,2,2,opt); %% mu= 1.16
K_dk=minreal(ss(K_dk)); %D-K Controller

L=minreal(Gp*K,0.000001);
Ln=minreal(G*K_dk,0.000001);
T=feedback(L,eye(2),-1);
T=minreal(T);%Complementary Sensitivity
S=feedback(eye(2),L);
S=minreal(S); %Sensitivity

% Nominal Stability - DK
N_dk=lft(P,K_dk);
N_dk=minreal(N_dk,0.000001);
eig_dk=max(real(eig(N_dk)));

% Nominal performance - DK
Nf_dk=frd(N_dk,omega);
blk=[2 4];
[mubnds,muinfo]=mussv(Nf_dk(5:8,5:6),blk,'c');
muNP_dk=mubnds(:,1);
[muNPinf_dk,muNPw_dk]=norm(muNP_dk,inf);

% Robust Stability - DK
blk=[1 1; 1 1;1 1;1 1];
[mubnds,muinfo]=mussv(Nf_dk(1:4,1:4),blk,'c');
muRS_dk=mubnds(:,1);
[muRSinf_dk,muRSw_dk]=norm(muRS_dk,inf);

%Robust Performance - DK
blk=[1 1; 1 1;1 1;1 1;2 4];
[mubnds,muinfo]=mussv(Nf_dk(1:8,1:6),blk,'c');
muRP_dk=mubnds(:,1);
[muRPinf_dk,muRPw_dk]=norm(muRP_dk,inf);

% Normalized SSV bounds for NP, RS, RP
figure(13)
loglog(muNP_dk,muRS_dk,muRP_dk);
title('\mu')
xlabel('Frequency (rad/s)');
legend('\mu_{NP}-DK','\mu_{RS}-DK','\mu_{RP}-DK');
ylim([0.07,2.5]);


% Time-Domain Simulations

% DK
L=minreal(G*K_dk,0.000001);
S=feedback(eye(2),L);
S=minreal(S);

Ad=A;
Bd=B(:,3); %Disturbance Channel
Cd=C;
Dd=0;
ltisys_d=ss(Ad,Bd,Cd,Dd);
Gd=minreal(zpk(ltisys_d));
H=feedback(G*K_dk,eye(2));

%Mixed-Sensitivity
Ln=minreal(G*K,0.000001);
S2=feedback(eye(2),Ln);
S2=minreal(S2);
H2=feedback(G*K,eye(2));

%Dist Rejection
figure(14)
step(minreal(S2*Gd));
hold on;
step(minreal(S*Gd));
title('Disturbance Rejection')
legend('Mixed','DK')
xlabel('Time');
xlim([0,170]);

%Reference Tracking
figure(15)
step(H2(:,1));
hold on;
step(H(:,1));
hold off;
title('Reference Tracking')
xlabel('Time');
xlim([0,170])
legend('Mixed','DK')

% Bode Plots

%Mixed-Sensitivity Design
figure(16)
mixed=G*K;
bodemag(mixed(:,1));

%DK
hold on;
dk=G*K_dk;
bodemag(dk(:,1));
legend('Mixed','DK')


% Sensitivity and Complementary Sensitivity

%Mixed-Sensitivity Design
loopmixed=loopsens(G,K);

%DK
loopdk=loopsens(G,K_dk);

%Sensitivity
figure(17)
smix=minreal(inv(eye(2)+G*K));
sdk=minreal(inv(eye(2)+G*K_dk));
bodemag(smix(:,1))
hold on;
bodemag(sdk(:,1))
legend('Mixed','DK')
title('Sensitivity');

%Complementary Sensitivity
figure(18)
csmix=minreal(G*K*(inv(eye(2)+G*K)));
csdk=minreal(G*K_dk*(inv(eye(2)+G*K_dk)));
bodemag(csdk(:,1))
hold on
bodemag(csmix(:,1));
hold off
legend('DK','Mixed')
title('Complementary Sensitivity');

