%% Part A - Linearize the mode
%Estimates from last time CL_min = 0.67, V_min = 3.3
clc 
clear all
V_min=3.3; CL_min=0.67;

syms Dose V_hat CL_hat t V_lin CL_lin

%d0=subs(C,[V,CL],[V_hat,CL_hat]);
%d1 = subs(diff(C,V),[V,CL],[V_hat,CL_hat]);

C_p = Dose/V_hat*exp(-CL_hat/V_hat*t);
pd_wrtV = diff(C_p, V_hat); %partial derivative with respect to beta_hat
pd_wrtCL = diff(C_p, CL_hat); %partial derivative with respect to alpha_hat
f0 = C_p -  pd_wrtV*V_hat - pd_wrtCL*CL_hat;
f1 = pd_wrtV; 
f2 = pd_wrtCL;
C_lin = f0+f1*V_lin+f2*CL_lin;%in the format Dr. K wants
disp('C_p linearized about the Estimates V_hat and CL_hat')
disp(C_lin)

%% Part B -- Estimate parameters about V_lin and C_lin
T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];


%now we subsistute the symbolics with the values
f0_new = subs(f0, Dose, 50);
f0_new = subs(f0_new, t, T_obs);
f0_new = subs(f0_new, CL_hat, CL_min);
f0_new = subs(f0_new, V_hat, V_min);
f0_sub = double(f0_new); %new fnew

%now we subsistute the symbolics with the values
f1new = subs(f1, Dose, 50);
f1new = subs(f1new, t, T_obs);
f1new = subs(f1new, CL_hat, CL_min);
f1new = subs(f1new, V_hat, V_min);
f1new = double(f1new);

%now we subsistute the symbolics with the values
f2new = subs(f2, Dose, 50);
f2new = subs(f2new, t, T_obs);
f2new = subs(f2new, CL_hat, CL_min);
f2new = subs(f2new, V_hat, V_min);
f2new = double(f2new);

%Solving equation using matrix form of normal equations
A = [sum(f1new.*f1new), sum(f1new.*f2new); sum(f2new.*f1new), sum(f2new.*f2new)];
b = [sum((C_obs - f0_sub).* f1new);sum((C_obs - f0_sub).*f2new)];
mat = vpa(inv(A)*b);
V_lin_ans = mat(1);
CL_lin_ans = mat(2);


disp('Estimate V_lin')
disp(V_lin_ans);
disp('Estimate CL_lin')
disp(CL_lin_ans);

%% Part C -- standard error
clinnew = subs(C_lin, Dose, 50);
clinnew = subs(clinnew, t, T_obs);
clinnew = subs(clinnew, V_hat, V_min);
clinnew = subs(clinnew, CL_hat, CL_min);
clinnew = subs(clinnew, CL_lin, CL_lin_ans);
clinnew = subs(clinnew, V_lin, V_lin_ans);
clinnew = double(clinnew);

OSS = sum(C_obs-clinnew).^2;
%non-weighted squares objective function
N = 9; %number of observations
p = 2; %number of parameters 
variance = OSS/(N-p);
cov_mat = variance*inv(A);
disp('Covariance Matrix about V_lin and CL_lin (Part E)') %covariance matrix for part E
disp(cov_mat)
se = sqrt(diag(cov_mat))';
SE_V_lin = se(1);
disp('Standard Error for V_lin');
disp(SE_V_lin);
SE_CL_lin = se(2);
disp('Standard Error for CL_lin')
disp(SE_CL_lin)


%% part D -- Calculate Confidence Intervals for V_lin and CL_lin
N = 9; %N
p = 2; %total number of parameters
DF = N-p; %degrees of freedom
ts  = tinv([0.025 0.975], DF); %intervals and degrees of freedom 
CI_V = V_lin_ans + ts*SE_V_lin;
disp('95% Confidence Interval for V_lin')
disp(CI_V)
CI_CL = CL_lin_ans + ts*SE_CL_lin;
disp('95% Confidence Interval for CL_lin')
disp(CI_CL)

%% part f -- estimate Kel


syms Dose V_lin CL_lin t V_hat CL_hat alpha_hat beta_hat alpha beta;
C_p = beta_hat*exp(-alpha_hat*t);
beta = Dose/V_lin;
alpha = CL_lin/V_lin;
pd_wrtbeta = diff(C_p, beta_hat); %partial derivative with respect to beta_hat
pd_wrtalpha = diff(C_p, alpha_hat); %partial derivative with respect to alpha_hat

%linearized model
g0 = C_p -  pd_wrtalpha*alpha_hat - pd_wrtbeta*beta_hat; %fit model
g1 = pd_wrtalpha; 
g2 = pd_wrtbeta;

C_lin = g0+g1*alpha+g2*beta; %in form that Dr. K asks for


T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24]; 
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];

g0_new = subs(g0, alpha_hat, CL_hat/V_hat); 
g0_new = subs(g0_new, beta_hat, Dose/V_hat);
g0_new = subs(g0_new, Dose, 50);
g0_new = subs(g0_new, t, T_obs);
g0_new = subs(g0_new, CL_hat, CL_min);
g0_new = subs(g0_new, V_hat, V_min);
g0_sub = double(g0_new); %new fnew

g1new = subs(g1, alpha_hat, CL_hat/V_hat);
g1new = subs(g1new, beta_hat, Dose/V_hat);
g1new = subs(g1new, Dose, 50);
g1new = subs(g1new, t, T_obs);
g1new = subs(g1new, CL_hat, CL_min);
g1new = subs(g1new, V_hat, V_min);
g1new = double(g1new);

g2new = subs(g2, alpha_hat, CL_hat/V_hat);
g2new = subs(g2new, beta_hat, Dose/V_hat);
g2new = subs(g2new, Dose, 50);
g2new = subs(g2new, t, T_obs);
g2new = subs(g2new, CL_hat, CL_min);
g2new = subs(g2new, V_hat, V_min);
g2new = double(g2new);


%Solving equation using matrix form of normal equations
A = [sum(g1new.*g1new), sum(g1new.*g2new); sum(g2new.*g1new), sum(g2new.*g2new)];
b = [sum((C_obs - g0_sub).* g1new);sum((C_obs - g0_sub).*g2new)];
mat = vpa(inv(A)*b);
alpha_lin = mat(1);
disp('Estimate of alpha (kel = CL/V)') %I named the secondary parameter alpha ...
disp(alpha_lin);
%since alpha = CL/V
beta_lin = mat(2);
%since beta = Dose/V
%we can solve for V beause we know beta and Dose


alpha_hat = alpha_lin;
beta_hat = beta_lin;
syms Dose V_lin alpha_hat beta_hat CL_lin
clinnew = subs(C_lin, Dose, 50);
clinnew = subs(clinnew, t, T_obs);
clinnew = subs(clinnew, V_lin, V_lin_ans);
clinnew = subs(clinnew, CL_lin, CL_lin_ans);
clinnew = subs(clinnew, alpha_hat, alpha_lin);
clinnew = subs(clinnew, beta_hat, beta_lin);
clinnew = double(clinnew);

OSS = sum(C_obs-clinnew).^2;
N = 9;
p = 2;
DF = N-p;
var2 = OSS/DF;
covar2 = var2*inv(A);

%% Part F -- Calculates the standard  error for Kel
syms CL V
a_hat = [V_lin_ans, CL_lin_ans]; %primary parameters that we already estimated
a = [V, CL]; %make vector that has symbolic parameters
g_alpha = CL/V; %this is my beta or g_alpha

cov_mat = cov_mat; %covariance matrix from my part C/E

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


z = [];
q = [];
for i=1:p;
    for j = 1:p;
        z = temp(i) * temp(j);
    end
    q = sum(z);
end

se = sqrt(diag(sum(q)*cov_mat))';
disp('Standard Error Estimate kel:')
disp(se)
