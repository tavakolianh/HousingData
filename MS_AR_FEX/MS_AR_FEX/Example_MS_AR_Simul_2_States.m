% Example Script for MS_AR_Sim.m (run it in the same directory - just press f5)

% This script will simulate a 2 states MS autoregressive model with some pre defined parameters
% You can change the coefficients for your own flavor

addpath('m_Files');

clear; 

nr=1000;        % Number of observations in simulation
advOpt.distrib='Normal';         % Distribution to use ('Normal' or 't' - default = 'Normal')

Coeff.p=[.8 .1 ; ...    % Transition matrix (this also defines the value of k)
         .2 .9 ];
         
Coeff.C(1,1)= .01;   % Setting up the mean constant at State 1
Coeff.C(1,2)=-.01;   % Setting up the mean constant at State 2

Coeff.AR(:,1)=[ .2   .1 -.1];   % Setting up the AR param at State 1 (you can increase/decrease it as wanted)
Coeff.AR(:,2)=[-.15 -.1  .1];   % Setting up the AR param at State 2 (you can increase/decrease it as wanted)   

Coeff.Std(1,1)=.02;  % Setting up the standard deviavion at State 1
Coeff.Std(1,2)=.01;  % Setting up the standard deviavion at State 2

ar=size(Coeff.AR,1);    % Number of lags in autoregressive component
k=length(Coeff.p);      % Number of states simulated (according to Coeff.p)

[Simul_Out]=MS_AR_Sim(nr,Coeff,k,advOpt.distrib); % Simulating the series

rmpath('m_Files');