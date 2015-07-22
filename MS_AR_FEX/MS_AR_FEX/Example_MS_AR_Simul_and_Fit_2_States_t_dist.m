% Example Script for MS_AR_Simul.m and MS_AR_Fit.m (run it in the same directory - just press f5)

% This script will first simulate a 2 states MS model with some pre defined parameters
% using MS_AR_Sim.m and then fit the model with MS_AR_Fit.m.

% You can change the coefficients for your own flavor

% You'll see that fmincon() is doing a pretty good job estimating the MS
% model (cheers for it!)

addpath('m_Files');

clear; 

nr=1000;                    % Number of observations in simulation
advOpt.distrib='Normal';    % Distribution to use ('Normal' or 't' - default = 'Normal')
advOpt.std_method='white';  % Method for standard error calculation ('white' or 'newey_west' - default = 'white')

Coeff.p=[.9 .2 ; ...    % Transition matrix (this also defines the value of k)
         .1 .8 ];       % Each collum of it should sum to 1
                        % Notes that the transition matrix output for fitting may 
                        % not be exactly equal to this one since the
                        % the fit function doesn't know which states is
                        % 1 2, it just tries to find them and
                        % separate. This means that the pii (prob os staying at state i) can be
                        % deslocated through the collums, but it
                        % should respect the structure at Coeff.p
         
Coeff.C(1,1)= .01;   % Setting up the mean constant at State 1
Coeff.C(1,2)=-.01;   % Setting up the mean constant at State 2

Coeff.v(1,1)= 3;   % Setting up the student distribution parameter at State 1
Coeff.v(1,2)= 3;   % Setting up the student distribution parameter at State 2

Coeff.AR(:,1)=[ .2  .1 -.1];   % Setting up the AR param at State 1 (you can increase/decrease it as wanted)
Coeff.AR(:,2)=[-.5 -.1  .1];   % Setting up the AR param at State 2 (you can increase/decrease it as wanted)   

Coeff.Std(1,1)=.02;  % Setting up the standard deviavion at State 1
Coeff.Std(1,2)=.05;  % Setting up the standard deviavion at State 2

ar=size(Coeff.AR,1);    % Number of lags in autoregressive component
k=length(Coeff.p);      % Number of states simulated (according to Coeff.p)

[Simul_Out]=MS_AR_Sim(nr,Coeff,k,advOpt.distrib); % Simulating the series

x=Simul_Out.Sim_x;  

figure(2);  % Making sure matlab plots both figures
[Spec_Output]=MS_AR_Fit(x,ar,k,advOpt);    % Fitting to get parameters

rmpath('m_Files');