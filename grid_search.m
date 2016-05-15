function [ xmin,ymin ] = grid_search( xlow, xup, N )
dx = (xup-xlow)/N;
xmin=xlow;
ymin=f(xmin);
for i = 1:N
    x = xlow + i * dx;
    y = f(x);
    if (y<ymin)
        xmin=x;
        ymin=y;
    end
end
    function[objective_function]=f(CL, V)
        objective_function = WSS(CL, V, T_obs, C_obs)
    end
end
