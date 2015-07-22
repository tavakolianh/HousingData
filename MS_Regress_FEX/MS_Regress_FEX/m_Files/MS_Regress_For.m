% Function for forecasting in t+1 an Markov Switching regression model
% estimated with MS_Regress_Fit.m
%
%   Input:  Spec_Output - Specification output from estimation (check
%   MS_Regress_Fit.m)
%           newIndepData - New information that has arrived (maybe lagged variables ?)
%
%   Output: meanFor - mean forecast at t+1
%           stdFor  - sigma forecast at t+1
%
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   PhD Student in finance ICMA/UK (Starting october 2007)
%   Created: June/2007

function [meanFor,stdFor]=MS_Regress_For(Spec_Out,newIndepData)

% Retrieving variables from Spec_Output

S=Spec_Out.S;
Coeff=Spec_Out.Coeff;
Std=Spec_Out.Coeff.Std;

filtProb=Spec_Out.filtProb;
k=Spec_Out.k;
p=Spec_Out.Coeff.p;

count_S=0;
count_nS=0;
for i=1:length(S)

    if S(i)==1
        count_S=count_S+1;
        indep_S(:,count_S)=newIndepData(:,i);
    else
        count_nS=count_nS+1;
        indep_nS(:,count_nS)=newIndepData(:,i);
    end
end

% BUG FIX
% Older versions didn't worked for when there was no non switching
% parameter). The next lines will fix that. 
% Thanks Mr. Panagiotis Papanastasiou-Ballis (15/01/2008)

if exist('indep_nS')==0 
    indep_nS=0;
end

for j=1:k
    meanFor_S(1,j)=indep_nS*Coeff.nS_param+indep_S*Coeff.S_param(:,j); % mean forecast for each state
    stdFor_S(1,j)=Std(j); % sigma forecast for each state
end

meanFor=meanFor_S*(p*filtProb(end,:)'); % mean forecast
stdFor=stdFor_S*(p*filtProb(end,:)');   % std forecast
