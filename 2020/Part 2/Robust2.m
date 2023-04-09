clc; clear all; close all; %VIVEK VARMA(5227828) and ZHEHAN LI(5217350)
%% ASSIGNMENT 2- Multi-Variable Mixed-Sensitivity

load('Assignment_Data_SC42145.mat') %Loading the Data
Ar=A;
Br=B(:, 1:2); %First 2 column of B matrix(beta, tau)
Cr=C; %(wr, z)
Dr=[0 0; 0 0];
[numr, denr]=ss2tf(Ar, Br, Cr, Dr, 1);
[numr2, denr2]=ss2tf(Ar, Br, Cr, Dr, 2);

sr1=zpk(tf(numr(1, :), denr)); %Input beta Output wr
sr2=zpk(tf(numr(2, :), denr)); %Input beta Output z
sr3=zpk(tf(numr2(1, :), denr2)); %Input tau Output wr
sr4=zpk(tf(numr2(2, :), denr2)); %Input tau Output z
sg=[sr1 sr3; sr2 sr4]; %Plant TFs
%% Question 1

G0=evalfr(sg,0*i); %w=0
RGA0=G0.*pinv(G0)';

G1=evalfr(sg,0.8*pi*i); %w=0.8*pi
RGA1=G1.*pinv(G1)';
%% Question 2

p2=pole(minreal((sg(1,1)*sg(2,2))-(sg(1,2)*sg(2,1)), 0.0001)); %Poles
z2=zero(minreal((sg(1,1)*sg(2,2))-(sg(1,2)*sg(2,1)), 0.0001)); %Zeros
%% Question 3

s=zpk('s');
Wp11=((s/1.8)+(0.8*pi))/(s+0.8*pi*10^-4);
%% Question 4

Wp=[Wp11 0; 0 0.2]; 
Wu=[0.01 0; 0 ((5*10^-3)*(s^2)+(7*10^-4)*(s)+(5*10^-5))/((s^2)+(14*10^-4)*(s)+(10^-6))];
%% Question 5

Pc=minreal(balreal(augw(sg, Wp, Wu))); %Generalized plant
sp=size(Pc);  %Number of states of Generalized Plant
%% Question 6+7

[k, cl, gam, info]=hinfsyn(Pc, 2, 2); %Design controller
k=minreal(k);
k=ss(k);
m=size(k.A, 1); %Number of states of the controller
ktf=tf(k);
ol=ss(minreal(ss(sg*k))); %Open-loop TF
L=(1+tf(ol(1,1)))*(1+tf(ol(2,2)))-(tf(ol(2,1))*tf(ol(1,2))); %det(I+ol)
figure(1)
nyquist(L); %Generlized Nyquist plot

pol=pole(ol); %Checking number of OL RHP poles P
zer=tzero(ol); %Checking Number of CL RHP poles Z=Number of OL RHP zeros
%% Question 8

figure(2)
step(feedback(ol, eye(2))) %Reference tracking step input response
x=feedback(ol, eye(2))
stepinfo(x(1,1))
dcgain(x(1,1))
stepinfo(x(2, 1))
dcgain(x(2,1))

Ad=A;
Bd=B(:,3); %Third column of B matrix(Disturbance)
Cd=C;
Dd=[0; 0];
gd=zpk(ss(Ad, Bd, Cd, Dd)); %Disturbance TF

y=zpk(feedback(eye(2), ol)*gd); %CLTF with Disturbance TF
y=minreal(y, 0.0001);

figure(3)
step(y) %Output Disturbance Rejection Step Response
stepinfo(y(1))
dcgain(y(1))
stepinfo(y(2))
dcgain(y(2))


%% 2.1 MIMO Weighting Design
A21=A;
B21=B(:, 1:3); 
C21=C(1,:); %First row of C matrix(wr)
D21=[0 0 0];
[num121, den121]=ss2tf(A21, B21, C21, D21, 1); % Beta control input
[num221, den221]=ss2tf(A21, B21, C21, D21, 2); % Tau control input
[num321, den321]=ss2tf(A21, B21, C21, D21, 3); % Wind speed input

sr121=zpk(tf(num121, den121));
sr221=zpk(tf(num221, den221));
sr321=zpk(tf(num321, den321));

sg21=[sr121 sr221 sr321];
Gd21=sg21(3); %Disturbance transfer function
sg21w=sg21(1:2);

%%Question 2
s=zpk('s');
Wp21=((s/1.5)+(0.1*2*pi))/(s+((0.1*2*pi)*10^-4)); %Error Weight with M=1.5, A=10^-4, %Freq=0.1Hz

%%Question 3
HPF=makeweight(0.01, 0.01*2*pi, 100); %High-Pass Filter with cross over at 0.01Hz
LPF=makeweight(50, 0.01*2*pi, 0.01); %Low-Pass Filter with cross over at 0.01Hz 
Wu21=[tf(HPF) 0; 0 tf(LPF)]; %Controller Weight

%%Question 1
Pc2=minreal(balreal(augw(sg21w, Wp21, Wu21))); %Generalized plant
%% Question 4
[k21, cl21, gam21, info21]=hinfsyn(Pc2, 1, 2); %Hinfinity controller
k21=minreal(k21);
k21=ss(k21);
m21=size(k21.A, 1);
ktf21=tf(k21);
ol21=ss(minreal(ss(sg21(1:2)*k21))); %Open loop system

figure(4) %Sensitivity Plot
bode(Gd21/(1+ol21));

figure(5) %Sensitivity vs Error Weight
bode(Gd21/(1+ol21));
hold on
bode(1/Wp21);
legend('Sensitivity Function','1/Wp')

figure(6) %Sensitivity
bode(ktf21*Gd21/(1+ol21));

figure(7) %Sensitivity vs Controller Weight
store=ktf21*Gd21/(1+ol21);
bode(store(1));
hold on;
bode(inv(Wu21(1, 1)));
legend('K*Sensitivity Function','1/Wu')

figure(8) %Sensitivity vs Controller Weight
bode(store(2));
hold on;
bode(inv(Wu21(2,2)));
legend('K*Sensitivity Function','1/Wu')

y21=zpk(feedback(eye(1), ol21)*sg21(3)); %Closed Loop system for disturbance rejection
y21=minreal(y21, 0.0001);
%% Question 5
% Refer to Simulink Model(trial.slx)
%This code needs to be run before running the Simulink model, since the
%transfer functions are imported into Simulink from this workspace.
%Rather than using the Scope to view, use the data inspector to view the
%outputs. They are more accurate. A manual switch has been placed so that
%we can switch between the Wind_Data provided and a noiseless sine wave
%easily. Double click on the switch to change the switch connection. Double
%click on the sine wave generator to change the parameters as desired. The
%variable names of the blocks are identical to the MATLAB code above.