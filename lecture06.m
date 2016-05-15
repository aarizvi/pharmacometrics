%% lecture 06

dsolve('Dy=2*y')

dsolve('Dy=2*y', 'y(0)=1')

syms A V Vmax Km
C = A/V

dsolve('DA = Vmax*C / Km+C')

%% one compartmental model 

syms ka kel dose F

[x,xa] = dsolve('Dxa=-ka*xa, Dx=ka*xa-kel*x', 'xa(0) = dose, x(0)=0')

