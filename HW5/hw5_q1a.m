%% Linearize equation
clc
clear all
%Estimates from last time
%V_min = 3.3, CL_min = 0.67
V_min = 3.3; CL_min=0.67;
syms Dose V_lin CL_lin t V_hat CL_hat alpha_hat beta_hat alpha beta;
C_p = beta_hat*exp(-alpha_hat*t);
beta = Dose/V_lin;
alpha = CL_lin/V_lin;
pd_wrtbeta = diff(C_p, beta_hat); %partial derivative with respect to beta_hat
pd_wrtalpha = diff(C_p, alpha_hat); %partial derivative with respect to alpha_hat

%linearized model
f0 = C_p -  pd_wrtalpha*alpha_hat - pd_wrtbeta*beta_hat;
f1 = pd_wrtalpha;
f2 = pd_wrtbeta;

C_lin = f0+f1*alpha+f2*beta;

%% Calculate estimates of V_lin and CL_lin
T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];

f0_new = subs(f0, alpha_hat, CL_hat/V_hat);
f0_new = subs(f0_new, beta_hat, Dose/V_hat);

f0_new = subs(f0_new, Dose, 50);
f0_new = subs(f0_new, t, T_obs);
f0_new = subs(f0_new, CL_hat, CL_min);
f0_new = subs(f0_new, V_hat, V_min);
f0_sub = double(f0_new); %new fnew

f1new = subs(f1, alpha_hat, CL_hat/V_hat);
f1new = subs(f1new, beta_hat, Dose/V_hat);
f1new = subs(f1new, Dose, 50);
f1new = subs(f1new, t, T_obs);
f1new = subs(f1new, CL_hat, CL_min);
f1new = subs(f1new, V_hat, V_min);
f1new = double(f1new);

f2new = subs(f2, alpha_hat, CL_hat/V_hat);
f2new = subs(f2new, beta_hat, Dose/V_hat);
f2new = subs(f2new, Dose, 50);
f2new = subs(f2new, t, T_obs);
f2new = subs(f2new, CL_hat, CL_min);
f2new = subs(f2new, V_hat, V_min);
f2new = double(f2new);

%Solving equation using matrix form of normal equations
A = [sum(f1new.*f1new), sum(f1new.*f2new); sum(f2new.*f1new), sum(f2new.*f2new)];
b = [sum((C_obs - f0_sub).* f1new);sum((C_obs - f0_sub).*f2new)];
mat = vpa(inv(A)*b);
alpha_lin = mat(1);
disp('Estimate of alpha (kel = CL/V)')
disp(alpha_lin);
%since alpha = CL/V
beta_lin = mat(2);
%since beta = Dose/V
%we can solve for V beause we know beta and Dose


%% Calculates standard error for estimates beta_hat and alpha_hat
alpha_hat = alpha_lin
beta_hat = beta_lin
V_lin_ans = V_lin;
CL_lin_ans = CL_lin;
syms Dose V_lin alpha_hat beta_hat CL_lin
clinnew = subs(C_lin, Dose, 50);
clinnew = subs(clinnew, t, T_obs);
clinnew = subs(clinnew, V_lin, V_lin_ans);
clinnew = subs(clinnew, CL_lin, CL_lin_ans);
clinnew = subs(clinnew, alpha_hat, alpha_lin);
clinnew = subs(clinnew, beta_hat, beta_lin);
clinnew = double(clinnew);

OSS = sum(C_obs-f0_sub-clinnew).^2;

%non-weighted squares objective function
N = 9; %number of observations
p = 2; %number of parameters 
variance = OSS/N-p;


cova = variance*(inv(A))
SE = sqrt(cova)

%% Calculates the confidence intervals for V_hat and CL_hat
syms CL V
a_hat = [V_lin_ans, CL_lin_ans]; %primary parameters that we already estimated
a = [V, CL]; %make vector that has symbolic parameters
g_alpha = CL/V; %this is my beta or g_alpha

covar = [0.3367, -0.0278; -0.0278, 0.0564] %this is just my covariance matrix from part C/E


p = 2; %number of parameters
temp = []; %storing the partial derivatives
pd = []; %partial derivatives
for i=1:p;
    pd = diff(g_alpha,a(i)); %partial derivatives
    for j=1:p;   
        pd = subs(pd, a(j), a_hat(j)); %substitute C_lin and V_lin into partial derivatives
    end
    temp(i) = pd; %collect substituted partial derivatives
end

x = temp(1)*temp(2)*covar;




%% Calculates the covariance matrix for V_hat and CL_hat

%% Calculates the estiamte of the secondary parameters k_el=CL/V, 
% and obtains the standrard error of that estimate 