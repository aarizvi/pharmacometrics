function [output] = rk4(F, x1_IC, x2_IC, step_size, tstart, tfinal)
%x1_ic = x1 initial condition
%x2_ic = x2 initial condition
%step_size = step size
%tstart = initial time
%tfinal = final time
%to run question 2 type in rk4[100, 90, 0.01, 0, 100)
F = str2func(F);
h = step_size;
t = tstart:h:tfinal;
u = zeros(2,numel(t));
u(:,1) = x1_IC; %x1 initial condition
u(:,2) = x2_IC; %x2 initial condition
t_ = []
for i = 2:numel(t)
    u_ = u(:,i-1); %initial conditions in 2x1 vector
    k1 = h*F(t_,u_);
    k2 = h*F(t_+0.5*h,u_+0.5*h.*k1);
    k3 = h*F(t_+0.5*h,u_+0.5*h.*k2);
    k4 = h*F(t_+h,u_+h.*k3);
    u(:,i) = u(:,i-1) + (k1+2*k2+2*k3+k4)/6;
    t_(:,i+1) = t(i-1) + h; %update time

end
output = [t;u]'
%plot(t,u)
%legend('Equation 1', 'Equation 2')
%ylim([-20 100])





