clear all;
close all
clc;
%HW1 Q4 - OVERLAYING TWO GRAPHS WITH GIVEN DATA IN SEMI-LOG
Tobs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
Cobs = [100, 61, 37, 22, 14, 8, 5, 3, 2, 1, 1];
Dose = 100; V=1; CL=0.5;
T = [0:0.1:10];


C = Dose/V*exp(-CL/V*T);


figure,semilogx(Tobs,Cobs,'ro',T,C,'b','LineWidth',2)
title('Semi-logarithmic graphs')
xlabel('Time')
ylabel('Concentration')
legend('Observations','Theoretical')
