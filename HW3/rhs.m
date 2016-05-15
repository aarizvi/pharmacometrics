function [ du ] = rhs( t,u )
%t = time span
%u = initial conditions
V=2; kel=0.3; kin=9; kout=0.1; Smax=10; SC50=10;
du = zeros(2,1);
du(1) = -kel*u(1);
du(2) = kin-kout*(1+Smax*(u(1)/V)/(SC50+u(1)/V))*u(2);
end

