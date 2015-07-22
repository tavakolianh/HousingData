% Example Script MS_Regress_Fit.m

% Script for simulating a 2 state MS regression model. Just press f5 to run
% it

addpath('m_Files'); % add 'm_Files' folder to the search path

clear; 

nr=1000;        % Number of observations in simulation
distrib='t';    % The distribution assumption

Coeff.p=[.8 .1 ; ...    % Transition matrix (this also defines the value of k)
         .2 .9 ];

Coeff.S=[1 0 0];  % Setting up which variables at indep will have switching effect

Coeff.v(1,1)=5;     % Setting up the distribution parameter for t distribution at state 1
Coeff.v(1,2)=10;    % Setting up the distribution parameter for t distribution at state 1

Coeff.nS_param(1,1)= .2;    % Setting up the coefficients at non switching parameters 
Coeff.nS_param(2,1)=-.2;    % Setting up the coefficients at non switching parameters 

Coeff.S_param(:,1)= .5;    % Setting up the coefficients at switching parameters 
Coeff.S_param(:,2)=-.3;    % Setting up the coefficients at switching parameters 

Coeff.Std(1,1)=.02;  % Setting up the standard deviavion of the model at State 1
Coeff.Std(1,2)=.01;  % Setting up the standard deviavion of the model at State 2

% The explanatory variables are going to be random normal ones, so the mean and std of
% them is needed. 

Coeff.indepMean=[.2 .5 .1];  % Setting up the mean of independent (explanatory) variables 
Coeff.indepStd= [.1 .2 .2];  % Setting up the mean of independent (explanatory) variables 

k=size(Coeff.p,1);  % getting the value of k, according to Coeff.p

[Simul_Out]=MS_Regress_Sim(nr,Coeff,k,distrib);

rmpath('m_Files');