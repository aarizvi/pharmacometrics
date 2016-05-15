%% Homework 4: Minimization and Elements of Statistics

% Write a MATLAB function that calculates the weighted sum of squares (See
% WSS.m)
% objective function 

%function[output] = WSS(CL, V, T_obs, C_obs)
T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];
V = 3.1020; %estimates for NCM (see nc_parameters_hw4.m)
CL = 0.6482; %estimates from NCM  (see nc_parameters_hw4.m)

%plot(T_obs, C_obs);
%lets estimate them using our NCM estimates
obj_function = WSS(CL, V, T_obs, C_obs);
disp('Weighted Sum of Squares using Non-Compartmental Estimates')
disp(obj_function)

%using our estimates from NCM we get a WSS of 0.1722


%% part b -- apply grid search methods to minimize WSS. Obtain estimates of V and CL for the following bounds:
% 0<V<10 L, accuracy of 0.1
% 0<CL<1 L/hr, accuracy of 0.01

%Estimates of the model parameter values are then determined by the
%absolute minimum of the objective function.
%for nonlinear models, minimization of the objective function must be done
%numerically -- NONLINEAR REGRESSION

%this one is for clearance
T_obs = [0, 1, 2, 4, 8, 12, 16, 20, 24];
C_obs = [15, 12.3, 10.1, 6.7, 3.0, 1.4, 0.6, 0.4, 0.1];
CL_UB = 1; V_UB=10; %upper bounds
CL_LB = 0; V_LB=0; %lower bounds
N=100; %number to divide intervals
dCL =(CL_UB-CL_LB)/N; %step size intervals for CL
dV = (V_UB-V_LB)/N; %step size intervals for dV
CL_min = CL_LB; V_min=V_LB; %minimum values to input into function
ymin=WSS(CL,V, T_obs, C_obs); %WSS with NCM as guesses .... we wanna minimize this
for i = 1:N
        for j = 1:N
            CL = CL_LB+i*dCL; %commence the search
            V = V_LB+j*dV;
            E = WSS(CL,V,T_obs,C_obs); 
        if(E<ymin)
            V_min = V;
            CL_min = CL;
            ymin=E;
        end
    end
end

disp('Minimum Volume')
disp(V_min)
disp('Minimum Clearance')
disp(CL_min)
disp('Weighted Sum of Squares with Minimums')
disp(ymin)



