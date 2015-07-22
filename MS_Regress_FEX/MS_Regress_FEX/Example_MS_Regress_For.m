% Example Script MS_Regress_For.m

% This script will define some variables in the .mat file and then use all
% of it, but the last row, to fit the model and use forecast at t+1 using the information in the
% the last row. 
% I know the independent aren't lagged variables. its just an example...

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path

load Example_FEX.mat    % load some Data. Change for your data here

dep=some_variables(1:end-1,1);                  % Defining dependent variable from .mat file
Constant=ones(length(dep),1);                   % Defining a constant vector in mean equation (just an example of how to do it)
indep=[Constant some_variables(1:end-1,3:4)];   % Defining some explanatory variables
k=3;                                            % Number of States
S=[1 0 0];                                      % Defining which ones from indep will have switching effect (in this case variable 1 (constant), only)
advOpt.distrib='Normal';                        % The distribution assumption ('Normal' or 't')
advOpt.std_method='white';              % Defining the method for calculation of standard errors ('white' or 'newey_west')(default='white')

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

newIndepData=[1 some_variables(end,3:4)]; % this is the new information that is feed to the forcasting function

[meanFor,stdFor]=MS_Regress_For(Spec_Out,newIndepData); % calling forecasting function

% printing results to screen

fprintf(1,['\nThe mean forecast at t+1 is ',num2str(meanFor)]);
fprintf(1,['\nThe sigma forecast at t+1 is ',num2str(stdFor),'\n']);

rmpath('m_Files');

