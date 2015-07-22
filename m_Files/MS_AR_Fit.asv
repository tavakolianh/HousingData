% Function for estimation of a Autoregressive Markov Switching model with k
% states and p lags (MS(k)-AR(p))
%
%   The models swicthes in the whole mean equation including AR coefficients and constant
%
%
%   Input:  x - Time series to be modelled
%           ar - Order of autoregressive component
%           k - Number of states
%           advOpt - A structure with advanced options. The field advOpt.distrib 
%                    is the distribution assumption ('Normal' or 't',
%                    default='Normal'). the field advOpt.std_method is the
%                    method for calculating the standard error of the
%                    fitted coefficients ('white' or robust standard
%                    errors by 'newey_west'). Default = 'white'.
%
%   Output: Spec_Output - A structure with followning fields:
%
%               LL - Log likelihood of fitted model
%
%               Probs - States probabilities over time (each collum represents
%               each state, ascending order).
%
%               Coeff - All estimated coefficients for each state.
%                      (AR parameters, standard deviation, transition
%                      matrix - each collum represents each state, ascending order).
%
%   The coefficient std is the model's standard deviation.
%
%   In the Spec_Output.Coeff.p, the row i, collum j, represents the
%   transition probability of being in state j in t-1 and changing to state
%   i on time t.
%
%   In the Coeff.Ar field, the row i, collum j, represents the AR coeff in lag i at
%   state j
%
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   Phd Student in finance ICMA/UK (Starting october 2007)
%   Created: June/2007

function [Spec_Output]=MS_AR_Fit(x,ar,k,advOpt)

[nr,nc]=size(x);

if nargin<3
    error('The function needs at least 3 arguments');
end

if nargin==3
    distrib='Normal';
    std_method='white';
else
    if isfield(advOpt,'distrib')==0
        distrib='Normal';
    else
        distrib=advOpt.distrib;
    end
    
    if isfield(advOpt,'std_method')==0
        std_method='white';
    else
        std_method=advOpt.std_method;
    end
end

if strcmp(distrib,'Normal')==0&&strcmp(distrib,'t')==0
    error('The distrib input should be ''Normal'' or ''t''');
end

if k==1
    error('If youre using k=1 (one state), you dont need this routine. Just use arma/garch toolbox');
end

if nc>1
    error('Are you trying to estimate a autoregressive model in a matrix ?? Does that make any sense ??');
end

if ar<1
    error('ar input should be higher than 1')
end

if nr<ar
    error('Your model doesnt make any sense, more lags than observations ??');
end

if nr<=20
    error('You need more than 20 observations in order to estimate the model.')
end

%   Setting up the param0 by estimating a simple Ar(p) model with garch toolbox.
%   If no toolbox is found, the function uses some ordinary guesses.

try
    
    spec=garchset('R',ar,'display','off');
    spec=garchfit(spec,x);

    param0_AR=[];
    param0_C=[];
    for i=0:k-1
        param0_AR=[spec.AR -param0_AR];  % Building Ar part for param0 (the next state is always with negative sign from the other)
        param0_C=[spec.C -param0_C];
    end

    switch distrib
        case 'Normal'
            param0=[repmat(std(x),1,k) param0_C param0_AR (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
        case 't'
            param0=[repmat(std(x),1,k) repmat(5,1,k) param0_C param0_AR (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
    end

catch

    spec.AR=repmat(.2,1,ar);
    spec.C=mean(x);

    param0_AR=[];
    param0_C=[];
    for i=0:k-1
        param0_AR=[spec.AR -param0_AR];  % Building Ar part for param0 (the next state is always with negative sign from the other)
        param0_C=[spec.C -param0_C];
    end
    
    switch distrib
        case 'Normal'
            param0=[repmat(std(x),1,k) param0_C param0_AR (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
        case 't'
            param0=[repmat(std(x),1,k) repmat(2,1,k) param0_C param0_AR (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
    end
end

% Lower and upper bounds at MS estimation

switch distrib
    case 'Normal'
        lB=[repmat(0,1,k)   repmat( -inf,1,k) repmat(-2,1,ar*k) repmat(0,1,k^2)];
        uB=[repmat(inf,1,k) repmat(  inf,1,k) repmat( 2,1,ar*k) repmat(1,1,k^2)];
    case 't'
        lB=[repmat(0,1,k)   repmat(0,1,k)   repmat( -inf,1,k) repmat(-2,1,ar*k) repmat(0,1,k^2)];
        uB=[repmat(inf,1,k) repmat(inf,1,k) repmat(  inf,1,k) repmat( 2,1,ar*k) repmat(1,1,k^2)];
end

warning('off');

global k_global;    % Im using this global here because I need to pass the value of k 
% to @confuneq and i haven't found other way around
% it, since such function should have only one input
% parameter

k_global=k;

% Estimation using constrained minimization with @confuneq

options=optimset('fmincon');  
options=optimset(options,'display','off');                    

[param]=fmincon(@(param)MS_AR_Lik(x,param,ar,k,distrib,1),param0,[],[],[],[],lB,uB,@confuneq,options);

[V]=getvarMatrix(x,param,ar,k,distrib,std_method);

param_std=sqrt(diag((V))); % the std vector is the square of the diagonal of V (var-cov matrix)

% After finding param, filter it to the data to get estimated output

[arg1,Spec_Output]=MS_AR_Lik(x,param,ar,k,distrib,0);

% Clearing global variable

clear global k_global;

% calculating the smothed probabilities

for i=1:nr
    Spec_Output.smoothProb(i,1:k)=(Spec_Output.Coeff.p*Spec_Output.filtProb(i,1:k)')';
end

% Plotting time varying probabilities

plot([Spec_Output.filtProb]);
xlabel('Time');
ylabel('Filtered States Probabilities');

for i=1:k
    States{i}=['State ',num2str(i)];
end

legend(States);

% Sending output to matlab's screen

Spec_Output.Coeff.Std_std(1,1:k)=param_std(1:k);

switch distrib
    case 't'
        Spec_Output.Coeff.v_std(1,1:k)=param_std(k+1:k+k);
        Spec_Output.Coeff.const_std(1,1:k)=param_std(k+k+1:k+2*k);
        for i=0:k-1
            Spec_Output.Coeff.AR_std(1:ar,i+1)=param_std(k+2*k+1+i*ar:k+k+k+ar+i*ar);
        end
    case 'Normal'
        Spec_Output.Coeff.const_std(1,1:k)=param_std(k+1:2*k);

        for i=0:k-1
            Spec_Output.Coeff.AR_std(1:ar,i+1)=param_std(2*k+1+i*ar:k+k+ar+i*ar);
        end

end

fprintf(1,'\n');
fprintf(1,'\n***** MS Optimizations terminated. *****\n\n');
fprintf(1,['Final log Likelihood: ',num2str(Spec_Output.LL),'\n']);
fprintf(1,['Number of parameters: ',num2str(Spec_Output.Number_Parameters),'\n']);
fprintf(1,['Distribution Assumption -> ',distrib,'\n']);
fprintf(1,['Method for standard error calculation -> ',std_method,'\n']);

fprintf(1,'\n-----> Final Parameters <-----\n\n');

for j=1:k
    fprintf(1,['Parameters in State ',num2str(j),':\n\n']);
    fprintf(1,['AR param   -> ',num2str(Spec_Output.Coeff.AR(:,j)'),'\n']);
    fprintf(1,['Std Errors -> ',num2str(Spec_Output.Coeff.AR_std(:,j)'),'\n']);
    fprintf(1,['Constant   -> ',num2str(Spec_Output.Coeff.const(1,j)'),'\n']);
    fprintf(1,['Std Errors -> ',num2str(Spec_Output.Coeff.const_std(1,j)'),'\n']);
    fprintf(1,['Std Dev    -> ',num2str(Spec_Output.Coeff.Std(1,j)'),'\n']);
    fprintf(1,['Std Errors -> ',num2str(Spec_Output.Coeff.Std_std(1,j)'),'\n']);

    switch distrib
        case 't'
            fprintf(1,['v param    -> ',num2str(Spec_Output.Coeff.v(1,j)'),'\n']);
            fprintf(1,['Std Errors -> ',num2str(Spec_Output.Coeff.v_std(1,j)'),'\n']);
        case 'Normal'
            fprintf(1,'\n');
    end

end

fprintf(1,'-----> Transition Probabilities Matrix <-----\n\n');

factor=1000;
formated_p=floor(Spec_Output.Coeff.p*factor)/factor;

formated_p(formated_p<0)=0;

disp(num2str(formated_p));
