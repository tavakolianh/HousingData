% Example Script for MS_AR_Fit_param0.m (run it in the same directory)
% this example will pass an personalized initial parameter vector for the
% fitting function. 

% OBS: So far, the routines were only adapted for use with NORMAL
% Distribution. If you set the t distribution, they will not work

clear;

addpath('m_Files');

load Example_Data.mat; % load .mat file

x=ret;                      % Time series from example.mat
advOpt.distrib='Normal';    % Distribution to use ('Normal' or 't' - default = 'Normal')
advOpt.std_method='white';  % Method for standard error calculation ('white' or 'newey_west' - default = 'white')

% Passing on your param0 to the function which will then format it for the
% fitting function

Spec.AR0(:,1)=[ .2  .3];  % param0 for Ar parameters at state 1
Spec.AR0(:,2)=[-.1 -.2];  % param0 for Ar parameters at state 2

Spec.const0(1,1)= .01;  % param0 for mean constant at state 1
Spec.const0(1,2)=-.02;  % param0 for mean constant at state 2

Spec.Std0(1,1)=.02;  % param0 for model's standard deviation at state 1
Spec.Std0(1,2)=.03;  % param0 for model's standard deviation at state 2

Spec.p0=[.8 .1; ... % param0 for the transition matrix
         .2 .9];
     
k=size(Spec.p0,1);      % get k
ar=size(Spec.AR0,1);    % get ar

param0=spec2param0(Spec,ar,k);  % Use spec2param0() for translating the Spec to param0 vector

[Spec_Output2]=MS_AR_Fit_param0_ver(x,ar,k,param0,advOpt); % Fit the model

rmpath('m_Files');