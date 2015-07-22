% Function for simulation of an Autoregressive Markov Switching model with
% 2 states and p lags (MS(k)-AR(p))
%
%   The models switches in the whole mean equation including AR coefficients and constant
%
%   Input:  nr      - Number of rows (number of time periods)
%           Coeff   - Coeff structure with all the coefficients
%                       Coeff.AR -> autoregressive components
%                       Coeff.Std-> standard deviation
%                       Coeff.C  -> constant in mean
%                       Coeff.p  -> Transition probabilities matrix
%
%                       In the coeff struct, each collum represent each
%                       state. Take a look at the Example_Script_Simul_2_States.m,
%                       just follow it and you'll find very easy to
%                       simulate your own flavor (value of coefficients).
%           k       - Number of States
%           dist    - The distribution to be used in estimation ('Normal'
%           or 't')
%
%   Output: Simul_Out - A structure with following fields:
%
%               Sim_x - Simulated time series
%
%               Coeff - Structure with all coefficients
%
%               States- The "true" simulated states through time
% %
%   Author: Marcelo Scherer Perlin
%   Email:  marceloperlin@gmail.com
%   Phd Student in finance ICMA/UK (Starting october 2007)
%   Created: August/2007
%
%   Fell free to use it and/or modify it for your own interest.
%
%   Any comments are welcome, if you have a suggestion that will significantly
%   improve the code, please write it and send it to me. If the changes are interesting,
%   I'll repost the file with full aknowledgements.

function [Simul_Out]=MS_AR_Sim(nr,Coeff,k,dist)

if nargin()==3
    dist='Normal';
end

if strcmp(dist,'Normal')==0&&strcmp(dist,'t')==0
    error('The dist input should be ''Normal'' or ''t''');
end

if nr<0
    error('Wow, a negative number of observations! Are you implying that you''ve traveled through time with a time machine and discovered the beggining of days? Can I borrow the time machine ?');
end

if any( (sum(Coeff.p)>1.0001)|(sum(Coeff.p)<.999) )
    error('The sum of each collum in Coeff.p should be equal to 1 (they are probabilities..)')
end

if any(Coeff.p<0+Coeff.p>1)
    error('The Coeff.p stores probabilities and they should be lower than 1 and higher than 0, unless you have a new theory about probabilities (it should be a very convincing one)');
end

if size(Coeff.AR,2)~=k
    error('The Coeff.AR has a different number of collum than the value of k.');
end

if size(Coeff.Std,2)~=k
    error('The Coeff.Std has a different number of collum than the value of k.');
end

if size(Coeff.C,2)~=k
    error('The Coeff.C has a different number of collum than the value of k.');
end

if any(Coeff.Std<0)
    error('The Coeff.Std is a standard deviation and should be positive. Ohhhhh wait, are you saying that you''ve created the most perfect model, where the error is replaced by a varying degree of certanty ? Thats nice, you should try publish it...');
end

switch dist
    case 't'
        if isfield(Coeff,'v')==0
            error('The t distribution requires the v parameter for each state. Check the examples scripts');
        end

        if size(Coeff.v,2)~=k
            error('The Coeff.v has a different number of collumns than the value of k.');
        end
end

ar=size(Coeff.AR,1);    % getting ar value (max number of autoregressive)
first_idx=ar+1;         % first idx for conditional mean construction


Rnd=rand(nr,1); % random seed for states transition (uniform probabilities)

% Prealocation of large matrices

Cond_mean=zeros(nr,k);
States=zeros(nr,k);

States(1:first_idx-1,1)=1;  % starting with state 1 for first obs

for i=first_idx:nr  % Loop for construction of states (maybe vectorize later ??)

    state_past=find(States(i-1,:)==1);

    if Rnd(i,1)<Coeff.p(state_past,state_past)

        States(i,state_past)=1; % when staying at last state

    else    % when changing to other states

        idx_other=find(States(i-1,:)==0);
        Prob2=Coeff.p(:,state_past);

        a=[Coeff.p(state_past,state_past) ; Prob2(idx_other)];

        cum_sum=cumsum(a);
        sorted=sort([cum_sum ; Rnd(i,1)]); % throw the prob at cumsum of other states to get
        % where it stands (where to switch)

        idx=find(Rnd(i,1)==sorted)-1;      % find index

        States(i,idx_other(idx))=1;        % change state

    end
end

for i=first_idx:nr % Loop for construction of conditional mean for each state (vectorize later using filter())

    Cond_mean_Past(1:i-1,1)=sum(Cond_mean(1:i-1,:).*States(1:i-1,:),2);   % full mean in the past

    for ik=0:k-1    % mean for each state

        switch dist
            case 'Normal'
                Cond_mean(i,ik+1)=Coeff.C(1,ik+1)+flipdim(Coeff.AR(:,ik+1),1)'*Cond_mean_Past(i-ar:i-1,1)+randn()*Coeff.Std(1,ik+1);
            case 't'
                Cond_mean(i,ik+1)=Coeff.C(1,ik+1)+flipdim(Coeff.AR(:,ik+1),1)'*Cond_mean_Past(i-ar:i-1,1)+trnd(Coeff.v(1,ik+1))*Coeff.Std(1,ik+1);
        end
    end
end

% Passing up to output structure

Sim_x=sum(Cond_mean.*States,2);

Simul_Out.Sim_x=Sim_x;
Simul_Out.Coeff=Coeff;
Simul_Out.States=States;

% Plotting simulated series

figure(1);
plot(Simul_Out.Sim_x);
legend('Simulated Autoregressive MS Series');
xlabel('Time');
ylabel('Simulated Autoregressive MS Series');