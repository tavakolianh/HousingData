% Function for calculating standard errors of MS_Regress_Fit
%
% The calculation (aproximation) of the first and second derivative of the likelihood
% function is done by central finite differences
%
% Two method were implemented here, the 'white' and also the 'newey_west'
% method. The vector output from the last one is sometimes called "robust standard errors".
%
% References:
%
% NEWEY, B., WEST, D (1987) ‘A Simple, Positive semidefinite, Heteroskedasticity
% and Autocorrelation Consistent Covariance Matrix’ Econometrica, ed. 55, p. 347-370
%
% WHITE, H. (1984) ‘Asymptotic Theory for Econometricians’ New York:
% Academic Press.

function [V]=getvarMatrix(dep,indep,param,k,S,distrib,method)

fprintf(1,'\n\nCalculating Standard Error Vector...');

% Definition of a small change for each parameter.
% This has been shamelessly copied from function hessian_2sided.m provided
% at Kevin Sheppard's GARCH toolbox (http://www.kevinsheppard.com/)

myDelta=eps.^(1/3)*abs(param);

n=length(dep); % number of parameters

disp=0; % no display at likelihood evaluation

% First Derivative calculation

for i=1:numel(param)

    m=param;
    m(i)=param(i)+myDelta(i);
    [sumlik1,Output,logLikVec1]=MS_Regress_Lik(dep,indep,param,k,S,distrib,disp);
    [sumlik2,Output,logLikVec2]=MS_Regress_Lik(dep,indep,m,k,S,distrib,disp);

    s(:,i)=(logLikVec2-logLikVec1)/(2*myDelta(i));

end

sum_s_Matrix=zeros(numel(param));
for idx=1:n
    s_Matrix{idx,1}=s(idx,:)'*s(idx,:);
    sum_s_Matrix=sum_s_Matrix+s_Matrix{idx,1};
end

OP_Matrix=sum_s_Matrix/n;   % outer product matrix

switch method

    case 'newey_west'

        s_sum_part=zeros(length(param));
        for j=1:n

            sum_der=zeros(length(param));
            for i=j+1:n

                der{i}=s(i,:)'*s(i-j,:);
                sum_der=sum_der+der{i};

            end

            myGamma=sum_der/n;

            s_sum_part=s_sum_part+(1-j/(n+1))*(myGamma+myGamma');

        end

        S_param=OP_Matrix+s_sum_part;   % Matrix for newey_west

end

% Hessian (second partial derivatives Matrix) Calculation

sde=zeros(numel(param));
for i=1:numel(param)

    for j=1:numel(param)

        m1=param;
        m2=param;

        m1(i)=m1(i)+myDelta(i);
        m2(j)=m1(j)+myDelta(j);

        [sumlik1]=MS_Regress_Lik(dep,indep,param,k,S,distrib,disp);
        [sumlik2]=MS_Regress_Lik(dep,indep,m1,k,S,distrib,disp);
        [sumlik3]=MS_Regress_Lik(dep,indep,m2,k,S,distrib,disp);

        % second derivative by central finite difference
        
        sde(i,j)=(sumlik2-2*sumlik1+sumlik3)/(myDelta(i)*myDelta(j));

    end
end

H=sde/n; % the Hessian

switch method
    case 'white'
        V=1/n*(inv(-H)*OP_Matrix*inv(-H));
    case 'newey_west'
        V=1/n*(inv(H)*S_param*inv(H));
end