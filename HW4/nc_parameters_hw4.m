%Given data
Dose = 50
T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];
y = C_obs; %y is C observations
x = T_obs; %x is time points recorded

%lambda_z
lambda_z = log(y(length(y)-1)/y(length(y)))/(x(length(x))-(x(length(x)-1)));
%Calculate AUC
results = 0;
for i = 1:(length(T_obs)-1) %the first six indexes only for non-residual AUC (no need to exptrapolate to 0 because 0 given)
    trap_rule = ((y(i+1) + y(i))*(x(i+1) - x(i)))/2;
    results(i) = trap_rule;
end
AUC = sum(results(:));

%log trapezoidal rule to compute residual AUC 
AUC_residual = y(length(y)) / lambda_z; %8 because that is last index of column vectors T and C_obs to infinity
AUC_tot = AUC + AUC_residual;
AUC_tot

%Calculate AUMC
results = 0;
for i = 1:(length(T_obs)-1) %the first six indexes only for non-residual AUMC (no need to exptrapolate to 0 because 0 given)
    trap_rule_m1 = ((x(i+1)*y(i+1) + x(i)*y(i)) * (x(i+1)-x(i))) / 2 ;
    results(i) = trap_rule_m1; %collect local AUCs in results vector
end
AUMC = sum(results(:)); %sum all local AUCs besides residual AUC
AUMC_residual = ((x(length(x))*y(length(y)))/lambda_z) + (y(length(y))/((lambda_z)^2)); %calculate residual using AUC trap rule
AUMC_tot = AUMC + AUMC_residual

%Estimate of V
%V_z = Dose/lambda_z*AUC
%V_ss = Dose*AUMC/AUC^2
V_z = Dose/(lambda_z*AUC_tot)
V_ss = (Dose*AUMC_tot)/(AUC_tot^2)
V = V_ss

%Estimate of CL
CL = Dose/AUC_tot

