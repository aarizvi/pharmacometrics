function [ dydt ] = odefun(t,y)
dydt = [-99*y(1)-49*y(2); 98*y(1)+48*y(2)];
end

