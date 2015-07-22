% Function for estimation of a general Markov Switching regression with k states
%
%   The big advantage in this function is that you can select which
%   variables will have MS effect and which won't with the use of input S.
%   With a few matrix manipulations you can include constants, autoregressive
%   parameters and any other non latent type of variable.
%
%   Input:  dep     - Dependent Variable (vector)
%           indep   - Independent variables (a.k.a explanatory variables),
%                     can be a vector (1 variable) or a matrix (many
%                     variables)
%           k       - Number of States
%           S       - This variable controls for where to include a Markov Switching effect.
%                     As example, suppose you have a indep matrix with 3 collums
%                     (3 explanatory variables). You don't want to have MS
%                     in all parameters, but just in variables 1 and 3. For that you just set
%                     S=[1 0 1] and the function will estimate the model according to you MS choice
%                     (value 1 for MS effect, value 0 for no MS effect).
%                     Easy isn't it?
%           advOpt  - A structure with advanced options. The field advOpt.distrib is the distribution assumption ('Normal' or 't',
%                     default='Normal'). the field advOpt.std_method is the
%                     method for calculating the standard error of the
%                     fitted coefficients ('white' or robust standard
%                     errors by 'newey_west'). Default = 'white'.
%
%   Output: Spec_Output - A structure with followning fields:
%
%               LL - Log likelihood of fitted model
%               Probs - States probabilities over time (each collumn represents
%               each state, ascending order).
%
%               Coeff - All estimated coefficients for each state.
%                      (regression parameters, standard deviation, transition
%                      matrix - each collumn represents each state, ascending order).
%
%   In the Spec_Output.Coeff.p, the row i, collumn j, represents the
%   transition probability of being in state j in t-1 and changing to state
%   i on time t.
%
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   PhD Student in finance ICMA/UK
%   Created: August/2007

function [Spec_Output]=MS_Regress_Fit(dep,indep,k,S,advOpt)

% Error checking lines

if nargin<4
    error('The function needs at least 4 arguments');
end

if nargin==4
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

if strcmp(std_method,'newey_west')==0&&strcmp(std_method,'white')==0
    error('The input advOpt.std_method should be either ''white'' or ''newey_west''');
end

if strcmp(distrib,'Normal')==0&&strcmp(distrib,'t')==0
    error('The distrib input should be ''Normal'' or ''t''');
end

if size(dep,2)>2
    error('The dep variable should be a vector, not a matrix');
end

if size(dep,1)~=size(indep,1)
    error('The number of rows at dep should be equal to the number of rows at indep');
end

if k<2||sum(S)==0
    error('k should be an integer higher than 1 and S shoud have at least a index with value 1. If you trying to do a simple regression just use regress() at stat toolbox') ;
end

if size(indep,2)~=size(S,2)
    error('The number of collums at indep should match the number of collums at S') ;
end

if sum((S==0)+(S==1))~=size(S,2)
    error('The S input should have only 1 and 0 values (those tell the function where to place markov switching effects)') ;
end

% Calculation of some primary variables (number of variables, how many of them
% switching etc..)

nr=size(dep,1);

n_indep=size(indep,2);
n_S=sum(S);
n_nS=n_indep-n_S;

count=0;
countS=0;

% Checking which parameters will have switching effect

for i=1:length(S)

    if S(i)==1
        countS=countS+1;
        S_S(countS)=i;
        indep_S(:,countS)=indep(:,i);
    else
        count=count+1;
        S_nS(count)=i;
        indep_nS(:,count)=indep(:,i);
    end
end

if n_nS~=0
    param0_indep_nS=regress(dep,indep_nS); % simple Ols for param0 of non switching variables
else
    param0_indep_nS=[];
end

param_ols_S=regress(dep,indep_S); % simple Ols for param0 of switcing variables

param0_indep_S=[];
for i=0:k-1
    param0_indep_S=[param_ols_S' -param0_indep_S]; % building param0 of switching variables (changing sign of coefficients)
end

% Building the whole param0, which will be then feeded to fmincon

switch distrib
    case 'Normal'
        param0=[repmat(std(dep),1,k) param0_indep_nS' param0_indep_S (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
    case 't'
        param0=[repmat(std(dep),1,k) repmat(5,1,k) param0_indep_nS' param0_indep_S (nonzeros(repmat(.1,k,k)+eye(k)*(1-k*.1)))'];
end

% Calculation of lower and upper bounds

switch distrib
    case 'Normal'
        lB=[repmat(0,1,k)   repmat(-inf,1,length([param0_indep_nS' param0_indep_S])) repmat(0,1,k^2)];
        uB=[repmat(inf,1,k) repmat( inf,1,length([param0_indep_nS' param0_indep_S])) repmat(1,1,k^2)];
    case 't'
        lB=[repmat(0,1,k)   repmat(0,1,k)   repmat(-inf,1,length([param0_indep_nS' param0_indep_S])) repmat(0,1,k^2)];
        uB=[repmat(inf,1,k) repmat(inf,1,k) repmat( inf,1,length([param0_indep_nS' param0_indep_S])) repmat(1,1,k^2)];
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

[param]=fmincon(@(param)MS_Regress_Lik(dep,indep,param,k,S,distrib,1),param0,[],[],[],[],lB,uB,@confuneq,options);

[V]=getvarMatrix(dep,indep,param,k,S,distrib,std_method);

param_std=sqrt(diag((V)));

% After finding param, filter it to the data to get estimated output

[sumlik,Spec_Output]=MS_Regress_Lik(dep,indep,param,k,S,distrib,0);

% Clearing global variable

clear global k_global;

% calculating smoothed probabilities

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

switch distrib
    case 'Normal'

        Spec_Output.Coeff.Std_std(1,1:k)=param_std(1:k);
        Spec_Output.Coeff.Std(1,1:k)=param(1:k);
        Spec_Output.Coeff.nS_param_std(1:n_nS,1)=param_std(k+1:k+n_nS); % This is the non switching parameters

        for i=0:k-1
            Spec_Output.Coeff.S_param_std(1:n_S,i+1)=param_std(1+k+n_nS+i*n_S:k+n_nS+n_S+n_S*i);    % This is the switching parameters
        end

    case 't'

        Spec_Output.Coeff.Std_std(1,1:k)=param_std(1:k);
        Spec_Output.Coeff.v_std(1,1:k)=param_std(k+1:k+k);
        Spec_Output.Coeff.nS_param_std(1:n_nS,1)=param_std(k+k+1:k+k+n_nS); % This is the non switching parameters

        for i=0:k-1
            Spec_Output.Coeff.S_param_std(1:n_S,i+1)=param_std(k+1+k+n_nS+i*n_S:k+k+n_nS+n_S+n_S*i);    % This is the switching parameters
        end

end

Spec_Output.Number_Parameters=length(param);
Spec_Output.advOpt.distrib=distrib;
Spec_Output.advOpt.std_method=std_method;

% Sending output to matlab's screen

fprintf(1,'\n\n***** MS Optimizations terminated. *****\n\n');
fprintf(1,['Final log Likelihood: ',num2str(Spec_Output.LL),'\n']);
fprintf(1,['Number of parameters: ',num2str(Spec_Output.Number_Parameters),'\n']);
fprintf(1,['Distribution Assumption -> ',distrib,'\n']);
fprintf(1,['Method for standard error calculation -> ',std_method,'\n']);


fprintf(1,'\n-----> Final Parameters <-----\n');

fprintf(1,'\n---> Non Switching Parameters <---\n');

if n_nS==0
    fprintf(1,'\nThere was no Non Switching Parameters. Skipping this result');
else
    for i=1:n_nS
        fprintf(1,['\n Non Switching Parameter at Indep collumn ', num2str(S_nS(i))]);
        fprintf(1,['\n      Value:     ', num2str(Spec_Output.Coeff.nS_param(i))]);
        fprintf(1,['\n      Std error: ', num2str(Spec_Output.Coeff.nS_param_std(i))]);

    end
end

fprintf(1,'\n\n--->   Switching Parameters   <---\n');

for j=1:k
    fprintf(1,['\n      State ', num2str(j)]);
    fprintf(1,['\n          Standard Deviation:     ', num2str(Spec_Output.Coeff.Std(1,j))]);
    fprintf(1,['\n          Std Error:              ', num2str(Spec_Output.Coeff.Std_std(1,j))]);

    switch distrib
        case 't'
            fprintf(1,['\n          Degrees of Freedom:      ', num2str(Spec_Output.Coeff.v(1,j))]);
            fprintf(1,['\n          Std Error:               ', num2str(Spec_Output.Coeff.v_std(1,j))]);
    end
end

for i=1:n_S
    fprintf(1,['\n\n Switching Parameters for Indep collumn ', num2str(S_S(i)),'\n']);

    for j=1:k
        fprintf(1,['\n      State ', num2str(j)]);
        fprintf(1,['\n          Value:     ', num2str(Spec_Output.Coeff.S_param(i,j))]);
        fprintf(1,['\n          Std error: ', num2str(Spec_Output.Coeff.S_param_std(i,j))]);
    end
end

fprintf(1,'\n\n---> Transition Probabilities Matrix <---\n\n');

factor=1000;
formated_p=floor(Spec_Output.Coeff.p*factor)/factor;

formated_p(formated_p<0)=0;

for i=1:k
    disp(['     ' num2str(formated_p(i,:))]);
end