%% part a
%to use rk4 ... rk4(F, x1, x2, h, tstart, tfinal)
%F = function name as string
%x1 = initial condition for parameter 1
%x2 = initial condition for parameter 2
%h = step size
%tstart = start time 
%tfinal = final time
%x(0) = 100, r(0) = 90, h=0.01, tstart=0, tfinal=100
[output] = rk4('genvec', 100, 90, 0.01, 0, 100)
t_rk4 = output(:,1);
x_rk4 = output(:,2);
r_rk4 = output(:,3);

plot(t_rk4,x_rk4)
hold on
plot(t_rk4,r_rk4)
legend('Equation 1', 'Equation 2')
ylim([-20 100])

%part b ... find solution using ode45 and compare to solution
[t u] = ode45(@genvec, [0 100], [100 90])
%call genvec, tspan = [0 100], initial conditions = [100 90]
figure
plot(t,u, t_rk4, x_rk4, t_rk4, r_rk4)
legend('Equation ode45X', 'Equation ode45R', 'Equation rk4X', 'Equation rk4R')
%many more step sizes in part a than part b....most likely represents a
%true estimate of the value

