%% Question 1
%part a
% Determine the matrix form for the following one comp. model w/
% first-order abs
% a bolus dose is injected to compartment X1 with bioavail F at time 0

syms ka kel dose F x1 x2
dx1 = -ka*x1;
dx2 = ka*x1 - kel*x2;
M = jacobian([dx1;dx2],[x1, x2]) %jacobian matrix


%% part b - find eigenvalues lambda 1, lambda 2 for matrix M
[V, lambda] = eig(M); %V = eigenvectors, lambda = eigenvalues

%grab lambda1 and lambda2 from lambda*I matrix
lambda1 = lambda(1); %eigenvalue 1
lambda2 = lambda(4); %eigenvalue 2

%grab eigenvectors (column vectors of V matrix)
v1 = V(:,1); %eigenvector 1 
v2 = V(:,2); %eigenvector 2

disp('lambda1 is:')
disp(lambda1)
disp('lambda2 is:')
disp(lambda2)
disp('eigenvector v1 is:')
disp(v1)
disp('eigenvector v2 is:')
disp(v2)
%% part c - verify that solution can be written in the vector form

x0 = [F*dose;0]; %describe beginning

%use cramers rule to solve for c1 and c2
c1 = det([x0, v2])/det(V);
c2 = det([v1, x0])/det(V);
syms t x1 x2
%now write it like how professor K shows in handout
x = c1*v1*exp(-lambda1*t) + c2*v2*exp(-lambda2*t)

% show that x is a soln to DE from part A
f = M*x;
int(f) %x is a solution to the DE from part A





%% part d
[x1, x2] = dsolve('Dx1=-ka*x1','Dx2=ka*x1-kel*x2','x1(0)=F*dose', 'x2(0)=0')
%solution is bateman function
%x1 and x2 are the same as the solutions for x

