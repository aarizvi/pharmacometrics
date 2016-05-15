function [ partial_auc ] = Partial_AUC( t1,t2 )
%Partial_AUC(t1,t2) calculates the partial AUC^{t2}_{t1} using the rule
%from question 3 for arbitrary: 0 \leq t2 < infinity

%The rule from question 3 uses the slope (k) from an exponential interpolation

if t1 < 0
    error('both inputs must be greater than or equal to 0')
    else if t2 < 0
        error('both inputs must be greater than or equal to 0')
    end
end

y1 = 866.67;
y2 = 1.07;
x1 = 0;
x2 = 12;

k =  (log(y1 / y2)) / (x2 - x1); %slope of exponential interpolation

syms x;
assume(x > 0);

f = @(x) y1*exp(-k*(x-x1)); %writing  exponential function with y1 and x1 

auc_t2_0 = integral(f, 0, t2); %integral from 0 < x < t2
auc_t1_0 = integral(f, 0, t1); %integral from 0 < x < t1

partial_auc = auc_t2_0 - auc_t1_0;  
end

