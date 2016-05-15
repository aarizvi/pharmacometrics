%%question 3
%Find the Taylor series expansion of the function f(x) = exp(x)
%x = 0
%Truncate the Taylor series at n = 1, 4, 8
%Plot the graphs of the truncated series vs x
%use the intervals -2 <= x <= 2 together with the graph f(x)
%what is the min value of n to approximate f(x) over this interval to get
%accuracy of 0.0001
clear all
clc
syms x
f = exp(x);
t1 = taylor(f, 'Order', 1); %truncate taylor series at n=1
t4 = taylor(f, 'Order', 4); %truncate taylor series at n=4
t8 = taylor(f, 'Order', 8); %truncate taylor series at n=8

ezplot(t1)
hold on
ezplot(t4)
ezplot(t8)
ezplot(f)

xlim([-2 2]);
ylim([-1 4]);

legend('Approximation of exp(x) up to 0(x^{1})',...,
    'Approximation of exp(x) up to 0(x^{4})',...,
    'Approximation of exp(x) up to 0(x^{8})',...,
    'exp(x)',...,
    'Location', 'Southeast');

title('Taylor Series Expansion of f(x) = exp(x)')
hold off


%approximate accuracy 
%if we take nth order approximation, error is
%\frac{f^{n+1}(c)}{(n+1)!}x^{n+1}
%where c is from [a,b]
%taking a = 0 and b = 1 this is less than
%  \frac{e}{n+1} \leq \frac{3}{(n+1!)}

3/factorial(8) %n = 7


