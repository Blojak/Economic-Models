close all;
clear;
clc;

% Paramter Calibration
% End Date ZLB
T =5;

par.r = [0.05, 0];       % ->5%
par.in = [0,0.05];
par.rho= 0.05;
par.kappa =1;

% exogenous variables

% Eigenvalues
A = [0, -1 ; 
     -par.kappa, par.rho];
EV = eig(A);
par.delta = EV(1,1); 
par.lambda = EV(2,1);

% Declare time index
t = 0:0.01:2*T;
t = t';
%% determination of ir:
% for t<T --> ir > 0
% for t>= --> ir = 0
N   = length(t);
T_index = find(t==T);
B = (1/(par.lambda-par.delta)) * [par.lambda^2 , - par.delta^2; par.lambda, -par.delta];

[SV1,SV2,SV3 ] = IRF( t, N ,T, par,B );

%% Figures
% Inflation
figure
plot(t, SV1(:,2)*100,'Linewidth',2)
hold on
plot(t, SV2(:,2)*100,'r','Linewidth',2)
hold on
plot(t, SV3(:,2)*100,'--r','Linewidth',1.5)
title('Inflation (%)')
xlim([0,10])
ylim([-20,10])
hold on
plot([0,10],[0,0],'k')
hold on
plot([T,T],[-20,10],'k')
print -depsc IRF_inflation;


% Output gap
figure
plot(t, SV1(:,1)*100,'Linewidth',2)
hold on
plot(t, SV2(:,1)*100,'r','Linewidth',2)
hold on
plot(t, SV3(:,1)*100,'--r','Linewidth',1.5)
title('Output Gap (%)')
xlim([0,10])
ylim([-20,10])
hold on
plot([0,10],[0,0],'k')
hold on
plot([T,T],[-20,10],'k')
print -depsc IRF_output;
