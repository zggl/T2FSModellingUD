function Results=IT2TSK_NeuroFM_RSVD_testExample1func(J_max,Iter)
DNum=121;
In=linspace(-10,10,DNum);
In_dim=size(In');
Out=zeros(In_dim);
n=In_dim(2);
for i=1:In_dim(1)
   Out(i,1)=sinc(In(i)/pi); 
end
N=DNum;
%surfc(X1,X2,Out1)

QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
%% 
Net(1).IterTrainMax=10;
Net(1).mm=0.5;
Net(1).nn=1-Net(1).mm;
Net(1).rho=0.01;            % thresh old
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
Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
Net(1).Gj=zeros(N,J_max);           % firing strength
%% UD parameters
coli=1;
s=n+1;
if s==2
    temp=minmax(In);
    min_ranges_p=temp(:,1)'*ones(1,2);  
    max_ranges_p=temp(:,2)'*ones(1,2); 
elseif s==3
    min_ranges_p=[-15,10,10]  
    max_ranges_p=[15, 30,20]  
end
[X_scaled,Xij]=UniformDesignWithScale(J_max,s,coli,min_ranges_p,max_ranges_p);
for ite=1:Iter
    Net=IT2TSKNeuroFM_RSVD(Net,QBFAStruc,X_scaled(:,1));
    Results(ite,:)=[Net(4).MSE,Net(4).NRMSE];
end
% Net(2).m{1}
% Net(2).sigma{1}
% Net(3).Cj