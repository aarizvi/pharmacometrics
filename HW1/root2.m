function [ output ] = root2( a, b, c )
%Solving quadratic equation numerical coefficients as input (a, b, c) 
syms x x1 x2

y = a*x^2 + b*x + c;
x = sqrt(b^2 - 4*a*c);
x1 = (-b + x) / (2*a)
x2 = (-b - x) / (2*a)
end

