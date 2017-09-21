function [Results,data_PS,N,Net]=IT2TSK_NeuroFM_RSVD_StockPricefunc(J_max,Iter,i)
% The stock prices refer to 10 aerospace companies and were daily collected from January 1988 to October 1991,
% which is taken from the Bilkent University Function Approximation Public Repository 
% \url{http://funapp.cs.bilkent.edu.tr/DataSets/}
clc
% J_max % the number of rules
% Iter % 
% i % the stock index
Dir='StockPriceData';
filename=strcat(Dir,'\StockPriceData',num2str(i),'.mat');
data1=load(filename,'datai');
[data,data_PS] = mapminmax(data1.datai,eps,1);
In=data(:,1:end-1);
In_dim=size(In);
n=In_dim(2);
Out=data(:,end);
N=In_dim(1);
%% 
Net(1).IterTrainMax=1;
Net(1).mm=0.5;
Net(1).nn=1-Net(1).mm;
Net(1).rho=0.05;            % threshold
Net(1).tau=0.01;            % threshold
Net(1).Data={In,Out};
Net(1).In_dim=In_dim;
Net(1).width = abs(max(In(:))-min(In(:)))/(4*J_max);        
Net(1).DataStat=[minmax(In);minmax(Out)];
Net(1).Numb=N;
Net(1).J_max=J_max;
Net(1).DataDim=[In_dim(2),size(Out,2)];
Net(1).Sj=zeros(1,J_max);           %  Size of cluster j
Net(1).Dim={n,[n,J_max],J_max,1};	% Dimension of network
Net(2).m={zeros(n,J_max),zeros(n,J_max)};           % mean vector UMF/LMF
Net(2).m_Len=n*J_max;
Net(2).sigma_Len=n*J_max;
Net(2).sigma0=0.3*ones(In_dim(2),1);
Net(2).sigma={ones(n,J_max),ones(n,J_max)};         % derivation vector
Net(3).Cj=zeros(J_max,In_dim(2)+1);           % height of cluster
Net(3).Cstar=zeros(1,J_max);
Net(3).BSVD.A=zeros(N,J_max); 
Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
Net(1).Gj=zeros(N,J_max);           % firing strength
QBFAStruc.bound=[min(In(:));max(In(:));[0.1;Net(1).width]];% bound of center m and width sigma
%% UD parameters
coli=1;
s=In_dim(2);
min_ranges_p=min(In);  
max_ranges_p=max(In); 
uns = UDall(J_max,s); % uniform design method
% save('uns.mat','uns'); 
% load('uns.mat','uns')
[X_scaled,Xij]=UDWithScale(J_max,s,coli,min_ranges_p,max_ranges_p,uns);
Net(4).DataNameId=1;
Net(4).mapmmPSFlag{Net(4).DataNameId}=0;
for ite=1:Iter
    Net=IT2TSKNeuroFM_RSVD(Net,QBFAStruc,X_scaled);
    Results(ite,:, i)=[Net(4).MSE,Net(4).NRMSE];
end
% Net(2).m{1}
% Net(2).sigma{1}
% sprintf('MSE: %4d, RMSE:  %4d, NRMSE: %4d ',Net(4).MSE,Net(4).RMSE,Net(4).NRMSE)
