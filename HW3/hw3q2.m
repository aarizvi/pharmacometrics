%% question 2 -- find all steady state solutions (baselines) of the following PK/PD model:

syms kel kin kout Smax SC50 X V R

%system of ODEs
dx = -kel*X
dr = kin-kout*(1+Smax*(X/V)/(SC50+X/V))*R

A = jacobian([dx;dr], [X,R]); %jacobian matrix

%steady state
X = solve(dx==0,X);
R = solve(dr==0,R);
A_ss = subs(A); %jacobian matrix at steady-state
eig(A_ss) %eigenvalues < 0 ... stable

%question 2b
V=2; kel=0.3; kin=9; kout=0.1; Smax=10; IC50=10;
B = subs(A); %substitute jacobian matrix A with new conditions
eig(B) %eigenvalues < 0 ... stable