clc;
clear all;
close all;
load('Assignment_Data_SC42145.mat') %loading the data
A1=A;
B1=B(:, 1); %First column of B matrix(beta)
C1=C(1, :); %First row of C matrix(wr)
%C1=-1*C(1, :); %First row of C matrix(wr)
D1=0;
[num1, den1]=ss2tf(A1, B1, C1, D1);
sys=tf(num1, den1);
s=zpk(sys)% plant
%% QUESTION 1:

%%plant zeros & poles
poles1=pole(s) %poles of plant
zeros1=zero(s) %zeros of plant
figure(1)
rlocus(s); %rootlocus of plant
figure(2)
pzmap(sys); %pole-zero map of plant

%bode plot of plant
figure(3)
margin(sys);
[gm, pm ,wcg, wcp]=margin(sys);

%%Sensitiviy, Complementary Sensitivity of plant
figure(4)
loops = loopsens(sys,1); 
bode(loops.Si,'r',loops.Ti,'b',loops.Li,'g');
legend('Sensitivity','Complementary Sensitivity','Open-loop tf');
title('Frequency Characteristics of the Plant')
%% QUESTION 2:
dcgain(feedback(sys,1)) % the steady state value of closed loop plant
stepinfo(feedback(sys,1)) %information about the step response
step(feedback(sys,1)) %step response plot

ms1=20*log10(getPeakGain(loops.Si)); %sensitivity peak value
mt1=20*log10(getPeakGain(loops.Ti)); %complementary sensitivity peak value

% Four P controllers
P=[-5, -10, -16,-17];
% Step responses of four P controllers
for i=1:length(P)
    p=P(i);
    figure(4+i)
    step(feedback(p*sys, 1));
    stepinfo(feedback(p*sys, 1))
end

% Bode plots of four P controllers
for i=1:length(P)
    figure(8+i);
    p=P(i);
    loops2=loopsens(p*sys, 1);
    margin(loops2.Li)
end

% Sensitivity  and Complementary Sensitivity plots of four P controllers
for i=1:length(P)
    figure(12+i);
    p=P(i);
    loops2=loopsens(p*sys, 1);
    bode(loops2.Si,'r',loops2.Ti,'b',loops2.Li,'g');
end

% Three I controllers
I=[-1.7388/s,-0.46311/s,-0.24956/s];
% Step responses of three I controllers
for i=1:length(I)
    p=I(i);
    figure(16+i)
    step(feedback(p*sys, 1));
    stepinfo(feedback(p*sys, 1))
end

% Bode plots of three I controller
for i=1:length(I)
    figure(19+i);
    p=I(i);
    loops2=loopsens(p*sys, 1);
    margin(loops2.Li)
end
% Sensitivity and Complementary Sensitivity plots of three I controllers
for i=1:length(I)
    figure(22+i);
    p=I(i);
    loops2=loopsens(p*sys, 1);
    bode(loops2.Si,'r',loops2.Ti,'b',loops2.Li,'g');
end

%% QUESTION 3:

% Time Constant Represenation of Plant transfer function
sys
s.DisplayFormat='time constant' 

%First PI controller
Ti=4.938;
s=zpk('s')
Kp=-1;
Ki=Kp/Ti;
pitf1=Kp+(Ki/s)

%Plots
dcgain(feedback(pitf1*sys,1)) 
stepinfo(feedback(pitf1*sys,1))
figure(26)
step(feedback(pitf1*sys,1))
figure(27)
margin(pitf1*sys)
figure(28)
loops2=loopsens(pitf1*sys, 1);
bode(loops2.Si,'r',loops2.Ti,'b',loops2.Li,'g');
loop=feedback(pitf1*sys,1);
[mag,phase,w]=bode(loop);
n=1;while 20*log10(mag(n))>-3;n=n+1;end
bandwidth1=w(n)

%Second PI controller
Ti=4;
s=zpk('s')
Kp=-1;
Ki=Kp/Ti;
pitf2=Kp+(Ki/s)

%Plots
dcgain(feedback(pitf2*sys,1)) 
stepinfo(feedback(pitf2*sys,1))
figure(29)
step(feedback(pitf2*sys,1))
figure(30)
margin(pitf2*sys)
figure(31)
loops2=loopsens(pitf2*sys, 1);
bode(loops2.Si,'r',loops2.Ti,'b',loops2.Li,'g');
loop=feedback(pitf2*sys,1);
[mag,phase,w]=bode(loop);
n=1;while 20*log10(mag(n))>-3;n=n+1;end
bandwidth2=w(n)

%Third(Final) PI controller
Ti=4;
s=zpk('s')
Kp=-1.1;
Ki=Kp/Ti;
pitf3=Kp+(Ki/s)

%Plots
dcgain(feedback(pitf3*sys,1)) 
stepinfo(feedback(pitf3*sys,1))
figure(32)
step(feedback(pitf3*sys,1))
figure(33)
margin(pitf3*sys)
figure(34)
loops2=loopsens(pitf3*sys, 1);
bode(loops2.Si,'r',loops2.Ti,'b',loops2.Li,'g');
loop=feedback(pitf3*sys,1);
[mag,phase,w]=bode(loop);
n=1;while 20*log10(mag(n))>-3;n=n+1;end
bandwidth3=w(n)
looptf=feedback(pitf3*sys,1)% final closed loop tf of the system
%% QUESTION 4:
% disturbance plant
A2=A;
B2=B(:,3); %Third column of B matrix)
C2=C1;
D2=0;
[num2, den2]=ss2tf(A2, B2, C2, D2);
sys2=tf(num2, den2);
s2=zpk(sys2)%tf of disturbance plant
zero(s2)
pole(s2)

% open-loop tf from previous question
Ti=4;
s=zpk('s')
Kp=-1.1;
Ki=Kp/Ti;
pitf=Kp+(Ki/s); %SISO PI controller
L=sys*pitf; %open-loop tf

% step response
E=inv(1+L)*sys2
figure(35)
step(E) %disturbance step input
stepinfo(E)
title('Step response of the system to output disturbance')
dcgain(E) 

%Proposing alternative controllers for checking patterns in step response
%for output disturbance rejection
Ti2=5;
s=zpk('s')
Kp=-1.1;
Ki2=Kp/Ti2;
pitf2=Kp+(Ki2/s); 
L2=sys*pitf2;
E2=inv(1+L2)*sys2
stepinfo(E2)
dcgain(E2)

Ti3=3;
s=zpk('s')
Kp=-1.1;
Ki3=Kp/Ti3;
pitf3=Kp+(Ki3/s); 
L3=sys*pitf3; 
E3=inv(1+L3)*sys2
stepinfo(E3)
dcgain(E3)

Ti=4;
s=zpk('s')
Kp4=-0.6;
Ki4=Kp4/Ti;
pitf4=Kp4+(Ki4/s); 
L4=sys*pitf4; 
E4=inv(1+L4)*sys2
stepinfo(E4)
dcgain(E4)

Ti=4;
s=zpk('s')
Kp5=-1.6;
Ki5=Kp5/Ti;
pitf5=Kp5+(Ki5/s);
L5=sys*pitf5;
E5=inv(1+L5)*sys2
stepinfo(E5)
dcgain(E5)

%Plotting the above step responses parallely
figure(36)
step(E); %K=-1.1, Ti=4
hold on;
step(E2); %K=-1.1, Ti=5
hold on;
step(E3); % K=-1.1, Ti=5
hold on;
step(E4); % K=-0.6, Ti=4
hold on;
step(E5); % K=-1.6, Ti=4
legend('K=-1.1, Ti=4', 'K=-1.1, Ti=5', 'K=-1.1, Ti=5', 'K=-0.6, Ti=4', 'K=-1.6, Ti=4');