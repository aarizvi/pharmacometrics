function dudt = genvec(t,u)
%t = time span
%u = initial conditions
V=3;kel=0.3;keo=0.2; 
Ap = u(1);
Ce = u(2);
Cp = Ap/V;
dudt = [-kel*Ap;
        keo*(Cp-Ce)];
end
