clear all
clc
syms V D t k_el;
V_hat = 1;
k_el_hat = 0.7;

C_lin = (D/V_hat)*exp(-k_el_hat*t) + (-D/V_hat^2)*exp(-k_el_hat*t)*(V-V_hat) + ((-D*t)/V_hat)*exp(-k_el_hat*t)*(k_el-k_el_hat);

[coefficients, terms] = coeffs(C_lin);
coefficients = double(coefficients);
terms = vpa(terms);
%looking at the form $C_lin = \alpha(t,D) + \beta(t,D)*V +
%\gamma(t,D)*k_{el}

alpha_coeffs = coefficients(3:4)
alpha_terms = terms(3:4) 

beta_coeff = coefficients(1)
beta_term = terms(1)


gamma_coeff = coefficients(2)
gamma_term = terms(2)

