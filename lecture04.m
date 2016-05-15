%%%% Lecture 04 - Elements of Calculus %%%%%%%%
%% Functions of one variable
x = 0:0.01:4;
y = exp(x);
plot(x,y)

%% Multivariable Functions
%z = x^2 + y^2 + 1
%for -1 <= x <= 1, -1 <= y <= 1
[X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);
Z = X.^2+Y.^2+1; %need dot operators for adding the exponents
surf(X,Y,Z)

%how to explore 2D and 3D plots


%% Limits -- made to define continuity


%% Derivatives
x = 2; h=0.0001;
(exp(-(x+h)) - exp(-x))/h 
%when you are in need to calculate derivatives numerically...(sensitivity analysis)

%% Rules of Differentiation
syms x a c k;
diff(c,x)

diff(x^a, x) %powerful in PK 

diff(exp(k*x), x) %kinetics

diff(log(x),x)

diff(log10(x),x)

%% way to evaluate your derivatives
syms x 
y = diff(exp(-x))
x=2
eval(y)
% there is a operator that will call your expression and take the value of your symbol and evaluate it
% allegedly its close to 'substitution'
%% Taylor Series Expansion
% derivative of nth order -- calculate derivative of derivative ...
% recursive
%interpretation of first derivative is tangent line
%interpretation of second derivative is concavity ... curvature
%second derivative = 0 ... inflection point
% some functions you can keep calculating infinimastly 

%% Tangent line
x=0:0.1:4;
y=exp(-x);
y_tangent=0.1353-0.1353*(x-2);
plot(x,[y;y_tangent])

%% Partial Derivatives

%in essence -- it is a derivative, but if put in context of multivariate
%functions -- you make it a function of one variable by fixing all other
%variables to the variable you want to differentiate and the other
%variables become a function of that variable.

















