function[output] = OSS(Emax, EC50, Ce)
Eobs = [1,10.4,95.6,101.6,97.9,94.8,97.6,94.9,81.0,69.6,40.1,24.2,7.4,6.6];
Epred = (Emax.*Ce)./(EC50+Ce);
obj_function = (Eobs-Epred).^2; %weighted squares objective function
output = sum(obj_function);
end
