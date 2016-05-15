%%%%%%%% PHC 504 %%%%%%%%%%
%%%%%% Final Project %%%%%%
%%%%%%% Abbas Rizvi %%%%%%%

%A bolus dose of 1000 ug of a drug was injected IV into a subject.
%The drug changes a pharmacodynamic biomarker E that was measured at the times listed in the table.
%The drug plasma concentration Cp was described by the one-compartment model:
%dAp/dt = -kel*Ap, where Cp = Ap/V
%with kel = 0.3/h and V = 3L
%The drug distributes to the effect compartment (biophase) according to the
%following equation:
% dCe/dt = keo*(Cp-Ce) with Ce(0) = 0, where keo = 0.2/h
%And Ce is the drug concentration in the biophase.
%The effect E is described by the Emax model:
% E = (Emax*Ce)/(EC50+Ce)
%Assume that the data are available as MATLAB column vectors T and Eobs
%% We will first solve the equations
clc 
clear all
syms kel V keo t Cp(t) Ce(t) Ap(t)
%Solve equations by 2x2 SOE
A = [-kel 0; (1/V)*keo -keo]; %matrix form equations
Y = [Ap; Ce];
B = [0;0];
eqn = diff(Y) == A*Y + B;
init_cond = Y(0) == [1000; 0]; %Ap(0) = 1000, Ce(0) = 0

[ApSol(t) CeSol(t)] = dsolve(eqn, init_cond);

%parameters that were given:
%kel = 0.3/h V = 3L keo=0.2
%Data given:
T = [0,1,2,3,4,5,10,15,20,25,30,35,40,45];
Eobs = [1,10.4,95.6,101.6,97.9,94.8,97.6,94.9,81.0,69.6,40.1,24.2,7.4,6.6];

ApSol = subs(ApSol, kel, 0.3);
ApSol = double(subs(ApSol, t, T));

CpSol = ApSol/V;
CpSol = double(subs(CpSol, V, 3));

CeSol = subs(CeSol, keo, 0.2);
CeSol = subs(CeSol, kel, 0.3);
CeSol = subs(CeSol, V, 3);
CeSol = double(subs(CeSol, t, T));
figure
plot(T, CeSol, T, CpSol)
legend('Ce', 'Cp')
%% Question 1 -- calculates AUEC over the observed time interval using linear trapezoidal rule
% lets look at plot
plot(T, Eobs);
ylabel('Effect, %');
xlabel('Time, h');
xlim([0 50]);

x = T;
y = Eobs;

%Calculate AUEC
results = 0;
for i = 1:(length(x)-1) 
    trap_rule = ((y(i+1) + y(i))*(x(i+1) - x(i)))/2;
    results(i) = trap_rule;
end
AUEC = sum(results(:));
disp('Area Under Effect Curve');
disp(AUEC);

%first lets calculate the slope
%lambda_z = log(y(length(y)-1)/y(length(y)))/(x(length(x))-(x(length(x)-1)));
% %log trapezoidal rule to compute residual AUC 
% AUEC_residual = y(length(y)) / lambda_z; %8 because that is last index of column vectors T and C_obs to infinity
% AUEC_tot = AUEC + AUEC_residual;
% AUEC_tot
%% Question 2
%Apply RK4 method to solve numerically the PK model equations and evaluates
%Cp(t) and Ce(t) at arbitrary time T
%genvec is a function that I wrote that contains the equations in general
%vector form
h = 0.01; %step size
t = 0:0.01:45; %time interval
u = zeros(2,numel(t));
u(1,1) = 1000; %Ap initial condition
u(2,1) = 0; %Ce initial condition

t_ = [];
for i = 2:numel(t)
    u_ = u(:,i-1); %initial conditions in 2x1 vector
    k1 = h*genvec(t_,u_);
    k2 = h*genvec(t_+0.5*h,u_+0.5*h.*k1);
    k3 = h*genvec(t_+0.5*h,u_+0.5*h.*k2);
    k4 = h*genvec(t_+h,u_+h.*k3);
    u(:,i) = u(:,i-1) + (k1+2*k2+2*k3+k4)/6;
    t_(:,i+1) = t(i-1) + h; %update time
end
output = [t;u]';
% question 3 -- Plot Cp vs t and Ce vs t curves for 0<t<45
t_rk4 = output(:,1);
cp_rk4 = output(:,2)./3;
ce_rk4 = output(:,3);
figure
plot(t_rk4, ce_rk4, t_rk4, cp_rk4)
legend('Ce rk4', 'Cp rk4')
xlabel('Time')
ylabel('Concentration')

% compare to ode45
[t p] = ode45(@genvec, [0 45], [1000 0]);
cp_ode45 = p(:,1)./3; %V = 3, Cp=Ap/V
ce_ode45 = p(:,2);

figure
plot(t, cp_ode45, t, ce_ode45, t_rk4, cp_rk4, t_rk4, ce_rk4)
xlim([0 45])
legend('Cp ode45', 'Ce ode45', 'Cp rk4', 'Ce rk4')
xlabel('Time')
ylabel('Concentration')
%rk4 and ode45 overlap
%% question 4
% Calculates the ordinary least squares objective function for the PD data
% OSS(Emax,EC50, CeSol) -- must enter in CeSol from above
% function OSS
OSS(100, 3, CeSol); %random guesses of parameter values .... but calculates an OSS for Emax model

%% question 5 - Minimizes OSS(Emax,EC50) using grid search
Emax_UB = 120; EC50_UB=120; %upper bounds
Emax_LB = 0; EC50_LB=0; %lower bounds
N = 1200; %number to divide intervals so its 0.1 accuracy
dEmax =(Emax_UB-Emax_LB)/N; %step size intervals for dEmax = 0.1 accuracy
dEC50 = (EC50_UB-EC50_LB)/N; %step size intervals for dEC50 = 0.1 accuracy
Emax_min = Emax_LB; EC50_min=EC50_LB; %minimum values to input into function
ymin=OSS(100, 50, CeSol); %initial guess
for i = 0:N
        for j = 0:N
            Emax = Emax_LB+i*dEmax; %commence the search
            EC50 = EC50_LB+j*dEC50;
            E = OSS(Emax,EC50, CeSol);
        if(E<ymin)
            Emax_min = Emax;
            EC50_min = EC50;
            ymin=E;
        end
        end
end

disp('Minimum Emax');
disp(Emax_min);
disp('Minimum EC50');
disp(EC50_min);
disp('Minimum Ordinary Least Squares Objective Function');
disp(ymin);
%% Question 6 -- Calculates the covariance and correlation matrices from the estimate Emax_hat, ec50_hat
% linearize about estimates
% use equations to get Emax, ec50

%grab Ce from answer above
Ce_ans = [0, 51.9417, 81.0056, 94.8280, 98.7565, 96.4995, 57.0321, 25.7854, 10.5579, 4.1232, 1.5702, 0.5896, 0.2195, 0.0814];
syms Emax_hat ec50_hat Ce emax_lin ec50_lin

E = (Emax_hat*Ce)./(ec50_hat+Ce);
pd_emax = diff(E, Emax_hat);
pd_ec50 = diff(E, ec50_hat);

f0 = E - pd_emax*Emax_hat - pd_ec50*ec50_hat;
f1 = pd_emax;
f2 = pd_ec50;

E_lin = f0 + f1*emax_lin + f2*ec50_lin;
%estimate parameters about emax_lin and ec50_lin
f0_new = subs(f0, Emax_hat, Emax_min);
f0_new = subs(f0_new, ec50_hat, EC50_min);
f0_new = double(subs(f0_new, Ce, CeSol));

f1_new = subs(f1, ec50_hat, EC50_min);
f1_new = double(subs(f1_new, Ce, CeSol));

f2_new = subs(f2, Emax_hat, Emax_min);
f2_new = subs(f2_new, ec50_hat, EC50_min);
f2_new = double(subs(f2_new, Ce, CeSol));

%Solving equation using matrix form of normal equations
A = [sum(f1_new.*f1_new), sum(f1_new.*f2_new); sum(f2_new.*f1_new), sum(f2_new.*f2_new)];
b = [sum((Eobs - f0_new).* f1_new);sum((Eobs - f0_new).*f2_new)];
mat = inv(A)*b;
Emax_lin_ans = double(mat(1));
EC50_lin_ans = double(mat(2));

disp('Estimate Emax_lin');
disp(Emax_lin_ans);
disp('Estimate EC50_lin');
disp(EC50_lin_ans);

% covariance matrix
E_linnew = subs(E_lin, Ce, CeSol);
E_linnew = subs(E_linnew, ec50_hat, EC50_min);
E_linnew = subs(E_linnew, Emax_hat, Emax_min);
E_linnew = subs(E_linnew, ec50_lin, EC50_lin_ans);
E_linnew = double(subs(E_linnew, emax_lin, Emax_lin_ans));

oss = OSS(Emax_lin_ans, EC50_lin_ans, CeSol);

N = 14; %number of observations
p = 2; %number of parameters
var = oss/(N-p);
cov_mat = var*inv(A);
disp('Covariance matrix');
disp(cov_mat);

%correlation matrix
varpar = diag(cov_mat); %variance of parameters
corr_mat = cov_mat./sqrt(prod(varpar)); 
disp('Correlation Matrix');
disp(corr_mat);

%Diaganol of correlation matrix is not 1 ... seems strange ... will try MATLAB corrcov 
disp('Correlation Matrix using built-in MATLAB function')
disp(corrcov(cov_mat));


%% Question 7 -- Calculate SE and CV% and 95% CIs for the estimates Emax_hat and EC50_hat
%Standard Error
se = sqrt(diag(cov_mat))';
SE_Emax_lin = se(1);
disp('Standard Error for Emax');
disp(SE_Emax_lin);
SE_ec50_lin = se(2);
disp('Standard Error for EC50');
disp(SE_ec50_lin);

%Coefficient of variation (relative S.E.)
%if S.E. are very small, CV can be high, tells % of error that contributes
%to estimate
CV_Emax = abs((SE_Emax_lin./Emax_lin_ans));
disp('Emax %Coefficient of Variation');
disp(CV_Emax * 100)

CV_EC50 = abs((SE_ec50_lin./EC50_lin_ans));
disp('EC50 %Coefficient of Variation');
disp(CV_EC50 * 100);

% 95% confidence intervals
N = 14; %N
p = 2; %total number of parameters
DF = N-p; %degrees of freedom
ts  = tinv([0.025 0.975], DF); %intervals and degrees of freedom 
CI_Emax = Emax_lin_ans + ts*SE_Emax_lin;
disp('95% Confidence Interval for Emax')
disp(CI_Emax)
CI_EC50 = EC50_lin_ans + ts*SE_ec50_lin;
disp('95% Confidence Interval for EC50')
disp(CI_EC50)
%% Question 8 -- Calculates goodness of fit metrics: OSS_min, r^2, AIC
Eobs = [1,10.4,95.6,101.6,97.9,94.8,97.6,94.9,81.0,69.6,40.1,24.2,7.4,6.6];
Ehat = (Emax_lin_ans.*CeSol)./(EC50_lin_ans+CeSol);
%OSS min
OSS_ans = OSS(Emax_lin_ans, EC50_lin_ans, CeSol);
disp('Least Ordinary Sum of Squares for Ehat');
disp(OSS_ans);

%Calculate R-squared
num = [];
denom1 = [];
denom2 = [];
for i = 1:length(Eobs)
    num(i) = (Eobs(i)-mean(Eobs))*(Ehat(i)-mean(Ehat));
    denom1(i) = (Eobs(i) - mean(Eobs)).^2;
    denom2(i) = (Ehat(i) - mean(Ehat)).^2;
end

r = sum(num)/sqrt(sum(denom1)*sum(denom2));
disp('R-squared of Eobs vs Ehat');
disp(r^2);

% Calculate AIC
% AIC = N * ln(WSSR) + 2*p
AIC = N*log(OSS_ans) + 2*p;
disp('AIC for Ehat');
disp(AIC);
%% Question 9 - Output diagnostic plots
Eobs = [1,10.4,95.6,101.6,97.9,94.8,97.6,94.9,81.0,69.6,40.1,24.2,7.4,6.6];
Ehat = (Emax_lin_ans.*CeSol)./(EC50_lin_ans+CeSol);

% Eobs vs t overlaid with Ehat vs. t
figure
plot(T, Eobs, T, Ehat)
legend('Eobs', 'Ehat')
xlabel('Time')
ylabel('Effect %')

%Plot Eobs vs Ehat
figure
plot(Eobs, Ehat)
xlabel('Eobs')
ylabel('Ehat')






