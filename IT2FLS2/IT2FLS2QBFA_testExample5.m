%% SECTION TITLE
clc,clear
%% praparing data

% In=linspace(0,1,DNum);
% In_dim=size(In');
% Out=zeros(In_dim);
% n=In_dim(2);
% for i=1:In_dim(1)
%    etat=0.2*rand;
%    Out(i,1)=0.9*In(i)/(1.2*In(i)^3+In(i)+0.3)+etat; 
% end

addpath(genpath('./.'));

DName={'AutoMPG',... % 1
       'bank',...% 2
       'diabetes',... % 3
       'Istanbul Stock Exchange',... % 4
       };
DataNameId      = 4;  
Net(4).DataNameId = DataNameId;
savename=strcat('IT2FLS2QBFA',DName{DataNameId},'.txt');
% 样本数为偶数 
% 1= 'sinc', 2='AutoMPG', 3='bank', 4= 'diabetes',5='triazines' 
% 6='NonLinearSysIdentify', 7= 'Sincxy'   
DataName = char(DName(DataNameId))
Iter=1E4;
RelAll=zeros(Iter,4);
for iter=1:Iter
    iter
switch DataName
    case  'AutoMPG'  % 1  392x8
        tic
        iter
        trainNum=300; 
        AutoMPG=load('auto-mpg.data')';
        Net(4).mapmmPSFlag{DataNameId}=1;
        rand_sequence=randperm(size(AutoMPG,2),trainNum);
        rand_seqtest=setdiff(1:size(AutoMPG,2),rand_sequence);
        [In,Net(4).traindata_PS{DataNameId}] = mapminmax(AutoMPG(1:end-1,rand_sequence),0,1);
        [Out,Net(4).trainlabel_PS{DataNameId}] = mapminmax(AutoMPG(end,rand_sequence),0,1);
        [testdata,Net(4).testdata_PS{DataNameId}] = mapminmax(AutoMPG(1:end-1,rand_seqtest),0,1);
        [testlabel,Net(4).testlabel_PS{DataNameId}] = mapminmax(AutoMPG(end,rand_seqtest),0,1); 
        J_max=11;
        QBFAStruc.J_max=J_max;
        QBFAStruc.EvalNum=1;             % Iter times of QBFA;
        In_dim=size(In');
        QBFAStruc.VarNumber=In_dim(2)*J_max+In_dim(2)*J_max; % size: m{1}+sigma{1}
        QBFAStruc.Nb=5;                  % The number of bacteria
        QBFAStruc.Ncs=5;                 % Number of chemotactic steps
        QBFAStruc.Ns=5;                  % Limits the length of a swim
        QBFAStruc.Nre=5;                 % The number of reproduction steps
        QBFAStruc.Ned=2;                  % The number of elimination-dispersal events
        QBFAStruc.Ped=0.2;               % The probabilty that each bacteria will be eliminated/dispersed
        QBFAStruc.IterMax=3.5E3;
        QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
        %% 
        Net(1).IterTrainMax=10;
        Net(1).mm=0.5;
        Net(1).nn=1-Net(1).mm;
        Net(1).rho=0.02;            % thresh old
        Net(1).tau=0.05;            % threshold
        Net(1).Data={In',Out'};
        Net(1).In_dim=In_dim;
        n=In_dim(2);
        N=In_dim(1);
        Net(1).width = 0.3;        
        Net(1).DataStat=[minmax(In);minmax(Out)];
        Net(1).Numb=trainNum;
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
        Net=IT2TSKNeuroFM_RSVD_FLS2(Net,QBFAStruc);
        [Net1]=PredictIT2_FLS2(size(testdata,2),testdata',testlabel',Net);                     
        RelAll(iter,:)=[Net(4).RMSE,Net1(4).RMSE,Net(4).NRMSE,Net1(4).NRMSE];% LS:Training (RMSE)   Testing (RMSE) 
        if iter==1
            iterone=1;
        end
        if ~exist(savename,'file')
            iterone=iter;
            dlmwrite(savename,RelAll(iter,:),'delimiter','\t','precision','%.4d','newline', 'pc')
        elseif mod(iter,50)==0
            itertwo=iter;
            dlmwrite(savename,RelAll(iterone+1:itertwo,:),'delimiter','\t','precision','%.4d','-append','newline', 'pc');
            iterone=itertwo;
        end      
        toc
    case  'bank' % 2
        tic
        bank=load('bank.data');  %8192x8
        bankIdent1=bank';
        trainNum=5E3;
        bankSize=size(bankIdent1,2);        
        rand_sequence=randperm(bankSize,trainNum);
        rand_seqtest=setdiff(1:bankSize,rand_sequence);
        Net(4).mapmmPSFlag{DataNameId}=1;
        [In,Net(4).traindata_PS{DataNameId}] = mapminmax(bankIdent1(1:end-1,rand_sequence),0,1);
        [Out,Net(4).trainlabel_PS{DataNameId}] = mapminmax(bankIdent1(end,rand_sequence),0,1);
        [testdata,Net(4).testdata_PS{DataNameId}] = mapminmax(bankIdent1(1:end-1,rand_seqtest),0,1);
        [testlabel,Net(4).testlabel_PS{DataNameId}] = mapminmax(bankIdent1(end,rand_seqtest),0,1);  
        J_max=7;
        QBFAStruc.J_max=J_max;
        QBFAStruc.EvalNum=1;             % Iter times of QBFA;
        In_dim=size(In');
        QBFAStruc.VarNumber=In_dim(2)*J_max+In_dim(2)*J_max; % size: m{1}+sigma{1}
        QBFAStruc.Nb=1;                  % The number of bacteria
        QBFAStruc.Ncs=1;                 % Number of chemotactic steps
        QBFAStruc.Ns=1;                  % Limits the length of a swim
        QBFAStruc.Nre=1;                 % The number of reproduction steps
        QBFAStruc.Ned=2;                  % The number of elimination-dispersal events
        QBFAStruc.Ped=0.2;               % The probabilty that each bacteria will be eliminated/dispersed
        QBFAStruc.IterMax=3.5E3;
        QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
        Net(1).IterTrainMax=10;
        Net(1).mm=0.5;
        Net(1).nn=1-Net(1).mm;
        Net(1).rho=0.02;            % thresh old
        Net(1).tau=0.05;            % threshold
        Net(1).Data={In',Out'};
        Net(1).In_dim=In_dim;
        n=In_dim(2);
        N=In_dim(1);
        Net(1).width = 0.3;        
        Net(1).DataStat=[minmax(In);minmax(Out)];
        Net(1).Numb=trainNum;
        Net(1).J_max=J_max;
        Net(1).DataDim=[In_dim(2),size(Out,2)];
        Net(1).Sj=zeros(1,J_max);           %  Size of cluster j
        Net(1).Dim={n,[n,J_max],J_max,1};	% Dimension of network
        Net(2).m={zeros(n,J_max),zeros(n,J_max)};           % mean vector UMF/LMF
        Net(2).m_Len=n*J_max;
        Net(2).sigma_Len=n*J_max;
        Net(2).sigma0=0.2*ones(In_dim(2),1);
        Net(2).sigma={ones(n,J_max),ones(n,J_max)};         % derivation vector
        Net(3).Cj=zeros(J_max,In_dim(2)+1);           % height of cluster
        Net(3).Cstar=zeros(1,J_max);
        Net(3).BSVD.A=zeros(N,J_max); 
        Net(1).Gj=zeros(N,J_max);           % firing strength
        Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
        Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
        Net=IT2TSKNeuroFM_RSVD_FLS2(Net,QBFAStruc);
        [Net1]=PredictIT2_FLS2(size(testdata,2),testdata',testlabel',Net);                     
        RelAll(iter,:)=[Net(4).RMSE,Net1(4).RMSE,Net(4).NRMSE,Net1(4).NRMSE];% LS:Training (RMSE)   Testing (RMSE)
        if iter==1
            iterone=1;
         end
         if ~exist(savename,'file')
            iterone=iter;
            dlmwrite(savename,RelAll(iter,:),'delimiter','\t','precision','%.4d','newline', 'pc')
        elseif mod(iter,50)==0
            itertwo=iter;
            dlmwrite(savename,RelAll(iterone+1:itertwo,:),'delimiter','\t','precision','%.4d','-append','newline', 'pc');
            iterone=itertwo;
        end
        iter
        toc
    case 'diabetes' % 3  
        % diabetes
        iter
        diabetes2_data;     %   randomly generate new training and testing data for every trial of simulation
        traindata1=load('diabetes_train')'; %576x8
        testdata1=load('diabetes_test')';   %192x8
        nn.mapmmPSFlag{DataNameId}=0;       
        In=traindata1(1:8,:);
        Out=traindata1(end,:);
        Net(4).mapmmPSFlag{DataNameId}=0;
        trainNum=size(Out,2);
        testdata=testdata1(1:8,:);
        testlabel=testdata1(end,:); 
        J_max=5;
        QBFAStruc.J_max=J_max;
        QBFAStruc.EvalNum=1;             % Iter times of QBFA;
        In_dim=size(In');
        QBFAStruc.VarNumber=In_dim(2)*J_max+In_dim(2)*J_max; % size: m{1}+sigma{1}
        QBFAStruc.Nb=5;                  % The number of bacteria
        QBFAStruc.Ncs=5;                 % Number of chemotactic steps
        QBFAStruc.Ns=5;                  % Limits the length of a swim
        QBFAStruc.Nre=5;                 % The number of reproduction steps
        QBFAStruc.Ned=2;                  % The number of elimination-dispersal events
        QBFAStruc.Ped=0.2;               % The probabilty that each bacteria will be eliminated/dispersed
        QBFAStruc.IterMax=3.5E3;
        QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
        Net(1).IterTrainMax=1;
        Net(1).mm=0.5;
        Net(1).nn=1-Net(1).mm;
        Net(1).rho=0.02;            % thresh old
        Net(1).tau=0.05;            % threshold
        Net(1).Data={In',Out'};
        Net(1).In_dim=In_dim;
        n=In_dim(2);
        N=In_dim(1);
        Net(1).width = 0.2;        
        Net(1).DataStat=[minmax(In);minmax(Out)];
        Net(1).Numb=trainNum;
        Net(1).J_max=J_max;
        Net(1).DataDim=[In_dim(2),size(Out,2)];
        Net(1).Sj=zeros(1,J_max);           %  Size of cluster j
        Net(1).Dim={n,[n,J_max],J_max,1};	% Dimension of network
        Net(2).m={zeros(n,J_max),zeros(n,J_max)};           % mean vector UMF/LMF
        Net(2).m_Len=n*J_max;
        Net(2).sigma_Len=n*J_max;
        Net(2).sigma0=0.2*ones(In_dim(2),1);
        Net(2).sigma={ones(n,J_max),ones(n,J_max)};         % derivation vector
        Net(3).Cj=zeros(J_max,In_dim(2)+1);           % height of cluster
        Net(3).Cstar=zeros(1,J_max);
        Net(3).BSVD.A=zeros(N,J_max); 
        Net(1).Gj=zeros(N,J_max);           % firing strength
        Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
        Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
        Net=IT2TSKNeuroFM_RSVD_FLS2(Net,QBFAStruc);
        [Net1]=PredictIT2_FLS2(size(testdata,2),testdata',testlabel',Net);                     
        RelAll(iter,:)=[Net(4).RMSE,Net1(4).RMSE,Net(4).NRMSE,Net1(4).NRMSE];% LS:Training (RMSE)   Testing (RMSE)
        if iter==1
            iterone=1;
        end
        if ~exist(savename,'file')
            iterone=iter;
            dlmwrite(savename,RelAll(iter,:),'delimiter','\t','precision','%.4d','newline', 'pc')
        elseif mod(iter,50)==0
            itertwo=iter;
            dlmwrite(savename,RelAll(iterone+1:itertwo,:),'delimiter','\t','precision','%.4d','-append','newline', 'pc');
            iterone=itertwo;
        end
        iter
        toc
    case 'Istanbul Stock Exchange' % 4:536x9
    % Attribute Information:
    % Stock exchange returns. Istanbul stock exchange national 100 index, Standard & poors 500 return index, 
    % Stock market return index of Germany, Stock market return index of UK, Stock market return index of Japan, 
    % Stock market return index of Brazil, MSCI European index, MSCI emerging markets index
        iter
        ISEdata
        trainNum=400;
        rand_sequence=randperm(size(NonLinSysIdent,2),trainNum);
        rand_seqtest=setdiff(1:size(NonLinSysIdent,2),rand_sequence);
        Net(4).mapmmPSFlag{DataNameId}=1;
        [In,Net(4).traindata_PS{DataNameId}] = mapminmax(NonLinSysIdent(2:end,rand_sequence),-1,1);
        [Out,Net(4).trainlabel_PS{DataNameId}] = mapminmax(NonLinSysIdent(1,rand_sequence),-1,1);
        [testdata,Net(4).testdata_PS{DataNameId}] = mapminmax(NonLinSysIdent(2:end,rand_seqtest),-1,1);
        [testlabel,Net(4).testlabel_PS{DataNameId}] = mapminmax(NonLinSysIdent(1,rand_seqtest),-1,1);
        J_max=3;
        QBFAStruc.J_max=J_max;
        QBFAStruc.EvalNum=1;             % Iter times of QBFA;
        In_dim=size(In');
        QBFAStruc.VarNumber=In_dim(2)*J_max+In_dim(2)*J_max; % size: m{1}+sigma{1}
        QBFAStruc.Nb=5;                  % The number of bacteria
        QBFAStruc.Ncs=5;                 % Number of chemotactic steps
        QBFAStruc.Ns=5;                  % Limits the length of a swim
        QBFAStruc.Nre=5;                 % The number of reproduction steps
        QBFAStruc.Ned=2;                  % The number of elimination-dispersal events
        QBFAStruc.Ped=0.2;               % The probabilty that each bacteria will be eliminated/dispersed
        QBFAStruc.IterMax=3.5E3;
        QBFAStruc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
        Net(1).IterTrainMax=1;
        Net(1).mm=0.5;
        Net(1).nn=1-Net(1).mm;
        Net(1).rho=0.02;            % thresh old
        Net(1).tau=0.05;            % threshold
        Net(1).Data={In',Out'};
        Net(1).In_dim=In_dim;
        n=In_dim(2);
        N=In_dim(1);
        Net(1).width = 0.4;        
        Net(1).DataStat=[minmax(In);minmax(Out)];
        Net(1).Numb=trainNum;
        Net(1).J_max=J_max;
        Net(1).DataDim=[In_dim(2),size(Out,2)];
        Net(1).Sj=zeros(1,J_max);           %  Size of cluster j
        Net(1).Dim={n,[n,J_max],J_max,1};	% Dimension of network
        Net(2).m={zeros(n,J_max),zeros(n,J_max)};           % mean vector UMF/LMF
        Net(2).m_Len=n*J_max;
        Net(2).sigma_Len=n*J_max;
        Net(2).sigma0=0.2*ones(In_dim(2),1);
        Net(2).sigma={ones(n,J_max),ones(n,J_max)};         % derivation vector
        Net(3).Cj=zeros(J_max,In_dim(2)+1);           % height of cluster
        Net(3).Cstar=zeros(1,J_max);
        Net(3).BSVD.A=zeros(N,J_max); 
        Net(1).Gj=zeros(N,J_max);           % firing strength
        Net(4).JbestVec=zeros(1,Net(1).IterTrainMax);
        Net(4).MSEVec=zeros(1,Net(1).IterTrainMax);
        Net=IT2TSKNeuroFM_RSVD_FLS2(Net,QBFAStruc);
        [Net1]=PredictIT2_FLS2(size(testdata,2),testdata',testlabel',Net);
        MCIterDataFLS2{1,iter}=Net1(4).ErrorSquareVec;
        RelAll(iter,:)=[Net(4).RMSE,Net1(4).RMSE,Net(4).NRMSE,Net1(4).NRMSE];% LS:Training (RMSE)   Testing (RMSE)
        if iter==1
            iterone=1;
        end
        if ~exist(savename,'file')
            iterone=iter;
            dlmwrite(savename,RelAll(iter,:),'delimiter','\t','precision','%.4d','newline', 'pc')
        elseif mod(iter,50)==0
            itertwo=iter;
            dlmwrite(savename,RelAll(iterone+1:itertwo,:),'delimiter','\t','precision','%.4d','-append','newline', 'pc');
            iterone=itertwo;
        end
        iter
    otherwise
        warning('Wrong datasets, exit.')
        exit
end
end
save('IT2FLS2ErrorSquareVec.mat','MCIterDataFLS2')


