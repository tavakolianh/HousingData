% Example Script for MS_AR_Fit.m (run it in the same directory)

clear;

addpath('m_Files');

load Example_Data.mat; % load .mat file

ar=4;                       % Number of lags in autoregressive component
k=2;                        % Number of states 
x=ret;                      % Time series from example.mat
advOpt.distrib='Normal';    % Distribution to use ('Normal' or 't' - default = 'Normal')
advOpt.std_method='white';  % Method for standard error calculation ('white' or 'newey_west' - default = 'white')

[Spec_Output]=MS_AR_Fit(x,ar,k,advOpt); % fit the model

rmpath('m_Files');