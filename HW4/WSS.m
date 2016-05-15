function[output] = WSS(CL,V,T_obs,C_obs)
Dose=50;
%C_p = Dose/V*exp(-CL/V*T_obs); %one-compartmental model
%W = 1./(C_obs.^2); %weights
C_p = Dose/V*exp(-CL/V*T_obs); %one-compartmental model
W = 1./((C_obs).^2); %weights
obj_function = W.*((C_obs-C_p).^2); %weighted squares objective function
output = sum(obj_function);
end
