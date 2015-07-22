close all
clear
adr=fileparts(pwd);
if ispc
    load ([adr '\ExcelFiles\HousingDataTotal.csv']);  % For Windows
else
    load ([adr '/ExcelFiles/HousingDataTotal.csv']);
end
logdistrict=log(HousingDataTotal);
[trenddistrict,cycdistrict]=hpfilter(logdistrict,450);
cycdistrict1=hpfilter(cycdistrict,1);

% Number of states
addpath('m_Files');


ar1=ones(21,1);
ar=1;% Number of lags in autoregressive component
k=2;
x=cycdistrict(:,1);              % Time series from example.mat
advOpt.distrib='Normal';         % Distribution to use ('Normal' or 't' - default = 'Normal')
advOpt.std_method='newey_west';  % Method for standard error calculation ('white' or 'newey_west' - default = 'white')

[Spec_Output]=MS_AR_Fit(x,ar,k,advOpt); % fit the model
rmpath('m_Files');
Spec_OutputTotal.(['Spec_Output' num2str(1)])=Spec_Output;
FiltProb1(:,1)=Spec_Output.filtProb(:,1);
for i=2:21
    [pacf,lags,bounds] = parcorr(cycdistrict(:,i),4);
    ar1(i)=max(lags(((pacf<bounds(2))+(pacf>bounds(1)))>0));
    
end

for i=2:21
    addpath('m_Files');
    if i==4
        ar1(i)=1;
    else
        ar=ar1(i);                            % Number of lags in autoregressive component
    end
    k=2;
    x=cycdistrict(:,i);              % Time series from example.mat
    advOpt.distrib='Normal';         % Distribution to use ('Normal' or 't' - default = 'Normal')
    advOpt.std_method='newey_west';  % Method for standard error calculation ('white' or 'newey_west' - default = 'white')
    
    [Spec_Output]=MS_AR_Fit(x,ar,k,advOpt); % fit the model
    rmpath('m_Files');
    FiltProb1(:,i)=Spec_Output.filtProb(:,1);
    Spec_OutputTotal.(['Spec_Output' num2str(i)])=Spec_Output;
end
FiltProb2=1-FiltProb1;
for i=2:21
    names{i-1} = ['District ' num2str(i-1)];
    figure(1)
    subplot(6,4,i-1)
    plot(FiltProb1(:,1), 'k', 'LineWidth', 1.5)
    hold on
    plot(FiltProb1(:,i), 'r', 'LineWidth', 1.5, 'MarkerSize', 2.5)
    title (names{i-1},'Fontsize',8 );
    hold off
    grid on
    axis([1 76 min(min(FiltProb1(:,i),FiltProb1(:,1))) max(max(FiltProb1(:,i),FiltProb1(:,1)))]);
    xticklabel_rotate([1:4:76], 90, {{'75:1'}; {'76:1'}; {'77:1'}; {'78:1'}; {'79:1'}; {'80:1'}; {'81:1'}; {'82:1'}; {'83:1'}; {'84:1'}; {'85:1'}; {'86:1'}; {'87:1'}; {'88:1'}; {'89:1'}; {'90:1'}; {'91:1'}; {'92:1'}; {'93:1'}},'Fontsize',9)
end
for i=1:76
    for j=1:21
        if FiltProb1(i,j)>0.5
            Dating(i,j)=1;
        else
            Dating(i,j)=0;
        end
    end
end

% figure(2)
% plot(cycdistrict1(:,1), 'k', 'LineWidth', 1.5)
% grid on
% axis([1 76 min((cycdistrict1(:,1))) max((cycdistrict1(:,1)))]);
% xticklabel_rotate([1:4:76], 90, {{'75:1'}; {'76:1'}; {'77:1'}; {'78:1'}; {'79:1'}; {'80:1'}; {'81:1'}; {'82:1'}; {'83:1'}; {'84:1'}; {'85:1'}; {'86:1'}; {'87:1'}; {'88:1'}; {'89:1'}; {'90:1'}; {'91:1'}; {'92:1'}; {'93:1'}},'Fontsize',9)
% for i=2:21
%     names{i-1} = ['District ' num2str(i-1)];
%     figure(3)
%     subplot(6,4,i-1)
%     plot(cycdistrict1(:,1), 'k', 'LineWidth', 1.5)
%     hold on
%     plot(cycdistrict1(:,i), 'r-o', 'LineWidth', 1.5, 'MarkerSize', 2.5)
%     title (names{i-1},'Fontsize',8 );
%     hold off
%     grid on
%     axis([1 76 min(min(cycdistrict1(:,i),cycdistrict1(:,1))) max(max(cycdistrict1(:,i),cycdistrict1(:,1)))]);
%     xticklabel_rotate([1:4:76], 90, {{'75:1'}; {'76:1'}; {'77:1'}; {'78:1'}; {'79:1'}; {'80:1'}; {'81:1'}; {'82:1'}; {'83:1'}; {'84:1'}; {'85:1'}; {'86:1'}; {'87:1'}; {'88:1'}; {'89:1'}; {'90:1'}; {'91:1'}; {'92:1'}; {'93:1'}},'Fontsize',9)
% end

% for i=1:21
%     if i<10
%         ii=[ '0' num2str(i)];
%     else
%         ii=num2str(i);
%     end
%     eval(['District' ii '=cycdistrict1(:,' ii ');'])
% end
clc
if ispc
    save([adr '\Outputs\Housing.mat']);
else
    save([adr '/Outputs/Housing.mat']);
end


