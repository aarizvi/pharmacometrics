function[objective_function]=f(CL,V)
    T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
    C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];
    objective_function = WSS(CL, V, T_obs, C_obs)
end
