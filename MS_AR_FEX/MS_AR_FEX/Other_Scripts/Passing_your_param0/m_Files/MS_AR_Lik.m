% Likelihood Function for MS(k) AR(p) Estimation

function [sumlik,Output,logLikVec]=MS_AR_Lik(x,param,ar,k,distrib,disp)

nr=length(x);

% Organizing Coeffs for each state
% The order in param, including each state, is: [std, constants, AR, Transition_Probs]

Coeff.Std(1,1:k)=param(1:k);

switch distrib
    case 't'
        Coeff.v=param(k+1:k+k);
        Coeff.const(1,1:k)=param(k+k+1:k+2*k);

        for i=0:k-1
            Coeff.AR(1:ar,i+1)=param(k+2*k+1+i*ar:k+k+k+ar+i*ar);
        end

        for i=0:k-1
            Coeff.p(i+1,1:k)=param(end-k^2+1+i*k:end-k^2+k+i*k);
        end

    case 'Normal'

        Coeff.const(1,1:k)=param(k+1:2*k);

        for i=0:k-1
            Coeff.AR(1:ar,i+1)=param(2*k+1+i*ar:k+k+ar+i*ar);
        end

        for i=0:k-1
            Coeff.p(i+1,1:k)=param(end-k^2+1+i*k:end-k^2+k+i*k);
        end
        
end

% Getting first observation in filter

first_idx=ar+1;

% Prealocation of large matrices

Cond_mean=zeros(nr,k);
E=zeros(nr,k);
n=zeros(nr,k);
f=zeros(nr,1);
e=zeros(nr,k);

% Setting up first obs for MS filter

E(1:first_idx-1,1:k)=1/k;

% Vectorized main engine

for j=1:k % Getting full matrices for each state
    AR(:,j)=filter(Coeff.AR(:,j),1,x); % Building AR part for each state using filter()
    Cond_mean(first_idx:end,j)=Coeff.const(1,j)+AR(first_idx-1:end-1,j); % Conditional mean for each state
    e(:,j)=x(:,1)-Cond_mean(:,j); % Error for each state

    switch distrib
        
        case 'Normal'
            n(:,j)=(1./(Coeff.Std(1,j)*sqrt(2*pi()))).*exp(-(e(:,j).^2)./(2*Coeff.Std(1,j)^2)); % Normal prob density at each state
        case 't'
            n(:,j)=( gamma(.5.*(Coeff.v(1,j)+1)) ) ./ ( (gamma(.5.*Coeff.v(1,j))).*sqrt(pi().*Coeff.v(1,j).*Coeff.Std(1,j).^2).* ...
                (1+(e(:,j).^2)./(Coeff.v(1,j).*Coeff.Std(1,j)^2)).^(.5.*(Coeff.v(1,j)+1)) );
    end

end

for i=first_idx:nr

    f(i,1)=ones(k,1)'*(Coeff.p*E(i-1,:)'.*n(i,:)'); % MS Filter equation
    E(i,:)=((Coeff.p*E(i-1,:)'.*n(i,:)')/f(i,1));   % MS Filter equation for probabilities

end

% Negative sum of log likelihood for fmincon (fmincon minimizes the function)

sumlik=-sum(log(sum(f(first_idx:end,:),2)));

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
Output.ar=ar;
Output.k=k;
Output.Number_Parameters=length(param);
Output.param=param;
Output.condMean=sum(Cond_mean.*E,2);
Output.condStd=sum(repmat(Coeff.Std,nr,1).*E,2);

if disp
    fprintf(1,['\nSum log likelihood for ' distrib ' distribution - MS(',num2str(k),')-Ar(',num2str(ar),')-->', num2str(-sumlik)]);
end