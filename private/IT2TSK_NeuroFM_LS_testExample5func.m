function [Results,Net1]=IT2TSK_NeuroFM_LS_testExample5func(J_max,Iter,In,Out,nn)

In_dim=size(In');
n=In_dim(2);
N=In_dim(1);
%surfc(X1,X2,Out1)

QBFAStruc.bound=[min(In(:));max(In(:));[0.3;1]];% bound of center m and width sigma
%% 
Net(1).IterTrainMax=1;
Net(1).mm=0.5;
Net(1).nn=1-Net(1).mm;
Net(1).rho=0.02;            % thresh old
Net(1).tau=0.05;            % threshold
Net(1).Data={In',Out};
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
Net(1).Gj=zeros(N,J_max);           % firing strength
Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
Net(4).DataNameId=nn.DataNameId;
Net(4).mapmmPSFlag{nn.DataNameId}=nn.mapmmPSFlag{nn.DataNameId};
if nn.mapmmPSFlag{nn.DataNameId}==1
    Net(4).traindata_PS{nn.DataNameId}=nn.traindata_PS{nn.DataNameId};
    Net(4).trainlabel_PS{nn.DataNameId}=nn.trainlabel_PS{nn.DataNameId};
end
%% UD parameters
coli=1;
s=n;

Net(1).temp=minmax(In);
min_ranges_p=Net(1).temp(:,1)';   
max_ranges_p=Net(1).temp(:,2)';   

[X_scaled,Xij]=UniformDesignWithScale(J_max,s,coli,min_ranges_p,max_ranges_p);
%% 
Net(1).N0=55;    % Online Sequence LS
Net(1).Block=1;
Results1=[];
for ite=1:Iter
    Net2{ite}=IT2TSKNeuroFM_LS(Net,QBFAStruc,X_scaled,Net(1).N0,Net(1).Block);
    Results1(ite,:)=[Net2{ite}(4).MSE,Net2{ite}(4).NRMSE];   
end
[Results,ind]=min(Results1);
Net1=Net2{ind};
% Net(2).m{1}
% Net(2).sigma{1}
% sprintf('MSE: %4d, RMSE:  %4d, NRMSE: %4d ',Net(4).MSE,Net(4).RMSE,Net(4).NRMSE)