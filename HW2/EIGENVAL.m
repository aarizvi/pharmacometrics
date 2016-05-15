function [ lambda ] = EIGENVAL( A )
% HW2 Question 2
%Write a function EIGENVAL(A) where A is a real 2x2 matrix that determines
%if there are two distinct eigenvalues of the matrix A, i.e. solves the
%following equation
if size(A,1) > 1
%"Characteristic Polynomial"
p = poly(A); %A is an n-by-n matrix, 
%returns the n+1 coefficients of the characteristic polynomial of the matrix, det(lambdaI-A)
lambda = roots(p);
else
    disp('Must be at least 2x2 matrix')
end

