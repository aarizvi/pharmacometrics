%% Midterm Project Rizvi, Due Wed March 9 1400h
%% A bolus dose of 1000mg was injected IV into a subject. Plasma concentrations were recorded at specific times
% Create vectors containing time points
clear all
close all
clc

T = [0 0.5 1 2 4 6 8 12];
C_obs = [866.67 382.23 173.10 167.35 58.79 21.58 6.11 1.07];
%% Question 1 - Plot Concentration vs. Time data in linear scales
figure
plot(T,C_obs, 'r', T,C_obs,'r*')
xlim([0 15])
title('Concentration vs. Time')
xlabel('Time')
ylabel('Concentration')
legend('Observations')
%% Question 2 - Calculate C_max and t_max

[C_max, idx] = max(C_obs);%C_max is the maximum concentration that drug achieves after admin.
C_max
t_max = T(idx); %Describes the time at which C_max is observed
t_max
%% Question 3- Calculates lambda_z based on the last two observed points

%find slope (lambda_z) using last two points
x1 = 8;
x2 = 12;
y1 = 6.11;
y2 = 1.07;

lambda_z = log(y1/y2)/(x2-x1); %lambda z is slope of exponential interpolation
lambda_z
%% Question 4 - Calculate AUC using the linear trapezoidal method
%trapezoidal rule -- AUC of y=f(x)is approximated by the area of a
%trapezoid

y = C_obs; %y is C observations
x = T; %x is time points recorded
results = 0;
for i = 1:6 %the first six indexes only for non-residual AUC (no need to exptrapolate to 0 because 0 given)
    trap_rule = ((y(i+1) + y(i))*(x(i+1) - x(i)))/2;
    results(i) = trap_rule;
end

AUC = sum(results(:));

%log trapezoidal rule to compute residual AUC 
AUC_residual = y(8) / lambda_z; %8 because that is last index of column vectors T and C_obs to infinity

AUC_tot = AUC + AUC_residual;
AUC_tot
%% Question 5 -- Calculates AUMC using the linear trapezoidal method 

N = 8; %last recorded data point is at 8th index of column vectors T and C_obs
results = 0;
for i = 1:6 %the first six indexes only for non-residual AUMC (no need to exptrapolate to 0 because 0 given)
    trap_rule_m1 = ((x(i+1)*y(i+1) + x(i)*y(i)) * (x(i+1)-x(i))) / 2 ;
    results(i) = trap_rule_m1; %collect local AUCs in results vector
end
AUMC = sum(results(:)); %sum all local AUCs besides residual AUC
AUMC_residual = ((T(8)*C_obs(8))/lambda_z) + (C_obs(8) / (lambda_z)^2); %calculate residual using AUC trap rule
AUMC_tot = AUMC + AUMC_residual
%% Question 6 -- Calcuates MAT and MRT

%Since elimination was linear ....
MAT = 0;
MAT

MRT = (AUMC / AUC_tot) - MAT;
MRT


%% Question 7 -- Calculates NCA parameters: CL, V_ss, t_half, V_z
dose = 1000; %dose is 1000mg

CL = dose / AUC;
CL

V_ss = (dose * AUMC_tot) / (AUC_tot)^2;
V_ss

t_half = log(2)/lambda_z;
t_half

V_z = dose/(lambda_z*AUC_tot);
V_z




