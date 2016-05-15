%% Differential Equations 2

k = 0.5; h=0.1;

%could also do t=[0:h:1]
%could also do C=t; C has same dimension as t
%C(1) is initial condition
C(1) = 100; t(1) = 0;
C(2) = C(1) - k*C(1) * h; t(2) = t(1) + h;
C(3) = C(2) - k*C(2) * h; t(3) = t(2) + h;
for i=1:10,
    C(i + 1) = C(i) - k*C(i)*h;
    t(i+1)=t(i)+h;
end

for i = 0:10,
    True_sol(i+1) = 100*exp(-k*i*h);
end

[t',C',True_sol']
plot(t,C)
%could also plot with exact solution
plot(t,C,t,100.*exp(-k.*t)) %true solution is red line
legend('euler', 'true solution')
%the approximation error accumulates proportinally to the step size


