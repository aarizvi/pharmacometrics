clear all
clc

syms y1 y2 y3 a
eqn1 = -(1+a)*y1 + y2-2*y3 == -3;
eqn2 = y1-(2+a)*y2 == -2;
eqn3 = 2*y1-(3+a)*y3 == -2;

[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [y1, y2, y3])

det(A)

p = [-1 -6 -14 -11];
roots(p) %these are the values of the parameter "a" that will have no unique solution in the systems of equations

