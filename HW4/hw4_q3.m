%%homework 4 - question 3
clear all
clc

r = normrnd(10,1,[1000,1]); %draws 1000 measurements from normal population N(0,1)

figure
histfit(r)
xlim([6 14]);
ylim([0 100]);
title('Histogram of Frequency Distribution from a N(10,1) Population')
mean(r);
std(r);

%plot histogram of frequency distribution for arithemetic means of 1000
%samples with N=4 drawn from a normal population N(10,1).

r1 = normrnd(10,1,[1000,1]);
r2 = normrnd(10,1,[1000,1]);
r3 = normrnd(10,1,[1000,1]);
r4 = normrnd(10,1,[1000,1]);

r_am = (r1 + r2 + r3 + r4)/4;
figure
histfit(r_am)
xlim([6 14]);
ylim([0 100]);
title('Frequency Distribution for Arithmetric Means of 1000 Samples (N=4)');
mean(r_am);
std(r_am);

%The arithmetic mean of the 1000 random variables ~ N(10,1) from 4
%populations is closer to the expected value of a normal distributon. This
%is as per the Law of Large Numbers theorem and Central Limit Theorem. We can see that this holds
%true as the standard deviation for the arithmetic mean population is much
%lower than that of the one sampling alone. 
