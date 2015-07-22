% Script for random param0
% It may take a while to finish in slow computers..

% If you get log likelihood as -Inf in the screen, that means that the curret param0 is just not
% good enough for convergence. Leave it running and eventually it will go to next simulation.

addpath('m_Files');
clear;

load Example_Data.mat; % Change here for you own data

n=5;    % Number of simulations (for values higher then 10 it may take a lot of time..)

ar=1;   % Value of ar
k=2;    % Number of states
x=ret;  % passing time series

LL=zeros(n,1);
for i=1:n
    
    rndparam0_Out=rndparam0(x,ar,k);   % Creates a random param0
    
    try 
        [Spec_Output{i}]=MS_AR_Fit_param0_ver(x,ar,k,rndparam0_Out);  % Fit using rnd param0 and save Spec
        LL(i,1)=Spec_Output{i}.LL; % Save log likelihood    
    catch
        LL(i,1)=-inf;   % cases where the fmincon gave an error
    end
    
end

disp(' ');
for i=1:length(LL)
    fprintf(1,['\nSimulation ',num2str(i),' --> Log likelihood = ',num2str(Spec_Output{i}.LL)]);
end
disp(' ');  

rmpath('m_Files');