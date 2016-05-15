function [ output ] = f( CL, V )

output = Dose/V*exp(-CL/V*T_obs)
end

