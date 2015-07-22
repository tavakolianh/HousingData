% Example Script MS_Regress_Fit.m

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path

load Example_FEX.mat    % load some Data. Change for your data here

dep=some_variables(:,1);                % Defining dependent variable from .mat file
Constant=ones(length(dep),1);           % Defining a constant vector in mean equation (just an example of how to do it)
indep=[Constant some_variables(:,3:4)]; % Defining some explanatory variables
k=3;                                    % Number of States
S=[1 0 0];                              % Defining which ones from indep will have switching effect (in this case variable 1 (constant), only)
advOpt.distrib='t';                     % The Distribution assumption ('Normal' or 't')
advOpt.std_method='white';              % Defining the method for calculation of standard errors ('white' or 'newey_west')(default='white')

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');