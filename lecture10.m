%Slide 9 code
scatter(LBM,Muscle_Strength)
xlabel('Lean Body Mass')
ylabel('Muscle Strength')
s = regstats(Muscle_Strength, LBM, 'linear', {'yhat', 'beta', 'covb', 'rsquare', 'mse', 'tstat'});
hold on
plot(LBM,s.yhat)
hold off

%regstat ... requires column with dependent variables (responses), and
%columns with explanatory variables (X -- lean body mass)
%you have to specify what model you're using ... 
%mse .... is s value that WK talks about
% t statistic will give you a pvalue and tell you if distn is diff than 0

disp(s.beta)
disp(s.covb)
sqrt(s.covb(1,1)) %standard error .. reporting standard error ... important because parameters are dependent
%S.e only matters if it relates to estimate
sqrt(s.covb(2,2)) %estimate ... 

%relative standard error
a = sqrt(s.covb(1,1)/1-s.beta(1)) %absolute value of s.beta(1)

%test for being different from 0 is useful ... confidence intervals
%interpretation of precision of estimates
%if they are not statistically different from 0 ... WK doesn't keep them
%FDA might accept a 10% estimate .. but any higher will probably be
%considered poor

disp(s.rsquare) %if you take square root of something < 1 it will get bigger
disp(s.mse)
disp(s.tstat.beta)
disp(s.tstat.se)
disp(s.tstat.t) %2 p-values ... intercept is insignificant .. which we see anyway based off of SE
%2nd pvalue ... reject hypothesis that slope is 0, with a very small pvalue
disp(s.tstat.pval)
disp(s.tstat.dfe)

%polyfit ... linear regression is 1st order polynomial ... might be nice.

