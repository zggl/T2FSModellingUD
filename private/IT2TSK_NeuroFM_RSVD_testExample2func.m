function Results=IT2TSK_NeuroFM_RSVD_testExample2func(J_max,Iter)
s=2;
n=s;
J_max=8;
DNum=21;
if DNum<J_max
    error('Trainning data size must be >= %d',J_max)
end
N=DNum^2;
%% praparing data
% 
x1=linspace(0,pi,DNum);
x2=linspace(0,pi,DNum);
[X1,X2]=meshgrid(x1,x2);
In=[X1(:),X2(:)];
In_dim=size(In);
Out1=X2.*sin(X1)+X1.*cos(X2);
%surfc(X1,X2,Out1)

QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
%% 
Net(1).IterTrainMax=10;
Net(1).mm=0.5;
Net(1).nn=1-Net(1).mm;
Net(1).rho=0.01;            % thresh old
Net(1).tau=0.05;            % threshold
Net(1).Data={In,Out1(:)};
Net(1).In_dim=In_dim;
Net(1).width = abs(max(In(:))-min(In(:)))/(4*J_max);        
Net(1).DataStat=[minmax(In(:)');minmax(Out1(:)')];
Net(1).Numb=N;
Net(1).J_max=J_max;
Net(1).DataDim=[In_dim(2),size(Out1,2)];
Net(1).Sj=zeros(1,J_max);           %  Size of cluster j
Net(1).Dim={s,[s,J_max],J_max,1};	% Dimension of network
Net(2).m={zeros(s,J_max),zeros(s,J_max)};           % mean vector UMF/LMF
Net(2).m_Len=s*J_max;
Net(2).sigma_Len=s*J_max;
Net(2).sigma0=0.3*ones(In_dim(2),1);
Net(2).sigma={ones(s,J_max),ones(s,J_max)};         % derivation vector
Net(3).Cj=zeros(J_max,In_dim(2)+1);           % height of cluster
Net(3).Cstar=zeros(1,J_max);
Net(3).BSVD.A=zeros(N,J_max); 

Net(1).Gj=zeros(N,J_max);           % firing strength
%% UD parameters
coli=1;
if s==2
    temp=minmax(In');
    min_ranges_p=temp(:,1)';  
    max_ranges_p=temp(:,2)';  
elseif s==3
    min_ranges_p=[-15,10,10]  
    max_ranges_p=[15, 30,20]  
end
[X_scaled,Xij]=UniformDesignWithScale(J_max,s,coli,min_ranges_p,max_ranges_p);
%% 
for ite=1:Iter
    Net=IT2TSKNeuroFM_RSVD(Net,QBFAStruc,X_scaled);
    Results(ite,:)=[Net(4).MSE,Net(4).NRMSE];
end
% Net=Predict(Net(1).DataDim(1),In,Out1(:),J_max,Net)
% Net(2).m{1}
% Net(2).sigma{1}
% Net(2).sigma{2}
% sprintf('%4d',Net(4).MSE)
% Net(4).RMSE
% sprintf('%4d',Net(4).NRMSE)