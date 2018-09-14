clear
clc

global lags psth

% load the PSTH curves
load('psth_curves.mat');

% choose the curves that you'd like to fit
lags = lags_normal;
psth = psth_pd;

% set the initial guess for the model parameters that you're looking for.
% The routine will then refine the search. Let as assume for now that the
% parameters are x = [ A gamma w d]
x0 = [10 0.9 0.3 7];

% search the parameter values that best fit the PSTH curve
x = fminsearch(@impulseresp_error,x0);

% build the LTI system with the identified parameters
A = x(1); gamma = x(2); w = x(3); d = round(x(4));
num = [0 A*gamma*sin(w) 0];
den = [1 -2*gamma*cos(w) gamma^2];
Ts  = 0.1 ;
Hsys = tf(num,den,Ts,'variable','z^-1','InputDelay',d);

% compute the impulse response of the transfer function
y = impulse(Hsys,lags);

% compute the error between the actual impulse response
err = norm(psth-y,2);

% plot the results
plot(lags,psth,'k',lags,y,'r');
fprintf('MSE: %f\n',err);
fprintf('params: A= %f\t gamma= %f\t w= %f\t d= %f\n',x);
