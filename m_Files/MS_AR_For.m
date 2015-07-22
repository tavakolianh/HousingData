% Function for forecasting in t+1 of an Autoregressive Markov Switching
% model estimated with MS_AR_Fit.m
%
%   Input:  Spec_Output - Specification output from estimation (check
%   MS_AR_Fit.m)
%           x - Time series to be modelled
%           distrib - distribution assumption (normal or t)
%
%   Output: meanFor - mean forecast at t+1
%           stdFor  - sigma forecast at t+1
%
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   Phd Student in finance ICMA/UK (Starting october 2007)
%   Created: June/2007
%
%   Fell free to use it and/or modify it for your own interest.
%
%   Any comments are welcome, if you have a suggestion that will significantly
%   improve the code, please write it and send it to me. If the changes are interesting,
%   I'll repost the file with full aknowledgements.

function [meanFor,stdFor]=MS_AR_For(Spec_Output,x,distrib)

%   Retrieving variables from Spec_Output
if nargin==2
    distrib='Normal';
end

param=Spec_Output.param;
p=Spec_Output.Coeff.p;
AR=Spec_Output.Coeff.AR;
Std=Spec_Output.Coeff.Std;
const=Spec_Output.Coeff.const;
k=Spec_Output.k;
ar=Spec_Output.ar;

% Passing up param vector to likelihood function (this fix the case
% where you estimated the model over past data and is just filtering to
% new data. This way the probabilities over time are recalculated

[sumlik,Output]=MS_AR_Lik(x,param,ar,k,distrib,0);

filtProb=Output.filtProb;

for j=1:k
    meanFor_S(1,j)=const(j)+flipdim(AR(:,j),1)'*x(end-ar+1:end); % mean forecast for each state
    stdFor_S(1,j)=Std(j); % sigma forecast for each state
end

meanFor=meanFor_S*(p*filtProb(end,:)'); % mean forecast
stdFor=stdFor_S*(p*filtProb(end,:)');   % std forecast
