% Likelihood Function for MS(k) Regression

function [sumlik,Output,logLikVec]=MS_Regress_Lik(dep,indep,param,k,S,distrib,disp_out)

%   Calculation of some preliminary variables

nr=length(dep);

n_indep=size(indep,2);
n_S=sum(S);
n_nS=n_indep-n_S;

count=0;
countS=0;

for i=1:length(S)

    if S(i)==1
        countS=countS+1;
        indep_S(:,countS)=indep(:,i);
    else
        count=count+1;
        indep_nS(:,count)=indep(:,i);
    end
end

if n_nS==0
    indep_nS=zeros(nr,1);
    Coeff.nS_param=0;
end

% Organizing Coeffs for each state

switch distrib
    case 'Normal'
        Coeff.Std(1,1:k)=param(1:k);
        Coeff.nS_param(1:n_nS,1)=param(k+1:k+n_nS); % This is the non switching parameters

        for i=0:k-1
            Coeff.S_param(1:n_S,i+1)=param(1+k+n_nS+i*n_S:k+n_nS+n_S+n_S*i);    % This is the switching parameters
        end

        for i=0:k-1
            Coeff.p(i+1,1:k)=param(end-k^2+1+i*k:end-k^2+k+i*k); % Transition matrix
        end

    case 't'

        Coeff.Std(1,1:k)=param(1:k);
        Coeff.v(1,1:k)=param(k+1:k+k);
        Coeff.nS_param(1:n_nS,1)=param(k+k+1:k+k+n_nS); % This is the non switching parameters

        for i=0:k-1
            Coeff.S_param(1:n_S,i+1)=param(k+1+k+n_nS+i*n_S:k+k+n_nS+n_S+n_S*i);    % This is the switching parameters
        end

        for i=0:k-1
            Coeff.p(i+1,1:k)=param(end-k^2+1+i*k:end-k^2+k+i*k); % Transition matrix
        end

end

Cond_mean=zeros(nr,k);
e=zeros(nr,k);
n=zeros(nr,k);

% Vectorized main engine

for i=0:k-1
    Cond_mean(:,i+1)=indep_nS*Coeff.nS_param+indep_S*Coeff.S_param(:,i+1); % Conditional Mean for each state
    e(:,i+1)=dep(:,1)-Cond_mean(:,i+1); % Error for each state

    switch distrib
        case 'Normal'
            n(:,i+1)=(1./(Coeff.Std(i+1)*sqrt(2*pi()))).*exp(-(e(:,i+1).^2)./(2*Coeff.Std(1,i+1)^2)); % Normal prob density at each state
        case 't'
            n(:,i+1)=( gamma(.5.*(Coeff.v(1,i+1)+1)) ) ./ ( (gamma(.5.*Coeff.v(1,i+1))).*sqrt(pi().*Coeff.v(1,i+1).*Coeff.Std(1,i+1).^2).* ...
                (1+(e(:,i+1).^2)./(Coeff.v(1,i+1).*Coeff.Std(1,i+1)^2)).^(.5.*(Coeff.v(1,i+1)+1)) );
    end
end

% Prealocation of large matrices

E=zeros(nr,k);
f=zeros(nr,1);

% Setting up first probs of E

E(1,1:k)=1/k;

for i=2:nr

    f(i,1)=ones(k,1)'*(Coeff.p*E(i-1,:)'.*n(i,:)'); % MS Filter equation
    E(i,:)=((Coeff.p*E(i-1,:)'.*n(i,:)')/f(i,1));   % MS Filter equation for probabilities

end

% Negative sum of log likelihood for fmincon (fmincon minimizes the function)

sumlik=-sum(log(f(2:end)));

% Control for nan, Inf, imaginary. This works much simillar to a slap
% in the face of fmincon when it outputs nan, inf or imaginary

if isnan(sumlik)||isreal(sumlik)==0||isinf(sumlik)
    sumlik=inf;
end

% Building Output structure

logLikVec=f;
Output.Coeff=Coeff;
Output.filtProb=E;
Output.LL=-sumlik;
Output.k=k;
Output.param=param;
Output.S=S;
Output.condMean=sum(Cond_mean.*E,2);
Output.condStd=sum(repmat(Coeff.Std,nr,1).*E,2);

if disp_out==1
    fprintf(1,['\nSum log likelihood for MS Regression -->', num2str(-sumlik)]);
end