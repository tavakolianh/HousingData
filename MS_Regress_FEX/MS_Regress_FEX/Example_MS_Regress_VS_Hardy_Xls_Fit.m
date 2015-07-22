% This script will compare the results from MS_Regress_Fit.m against the
% results from the markov switching excel spreadsheet kindly provided by Professor Mary
% Hardy and Peter Hardy. You can download it at:
%
% http://www.soa.org/research/other-research-projects/tools-and-software/research-regime-switching-equity-model-workbook-version-1-0-updated-8-3-2004.aspx
%
% Be warned that the vba routine will not work in all versions of excel software.
%
% The model from the xls file is a markov swithcing with one constant at mean equation
% and also switching at volatitility with 2 states. That is:
%
% y(t)=m(st,t)+s(st,t)*e(t)
% e(t)~N(0,1)
%
% With a markov chain as 
% [p11 p21
%  p12 p22]
%
% The results from their spreadsheet is:
% m1	0.012540995     
% m2	-0.021898115
% s1	0.035565703
% s2	0.077259815
% p12	0.03449453
% p21	0.395653344

clear;

addpath('m_Files'); % add 'm_Files' folder to the search path

load Hardy_xls_data.mat    % load professor Hardy's Data

dep=log(Hardy_xls_Data(2:end)./Hardy_xls_Data(1:end-1));    % Changing prices to log returns
Constant=ones(length(dep),1);           % Defining a constant vector in mean equation (just an example of how to do it)
indep=Constant;                         % Defining a constant vector
k=2;                                    % 2 states in the model
S=1;                                    % Switching in the constant
advOpt.distrib='Normal';                % Normal distribution
advOpt.std_method='white';              % Defining the method for calculation of standard errors ('white' or 'newey_west')(default='white')

[Spec_Out]=MS_Regress_Fit(dep,indep,k,S,advOpt); % Estimating the model

rmpath('m_Files');