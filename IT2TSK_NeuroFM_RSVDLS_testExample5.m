clc
clear all
addpath(genpath('./.'));

DName={'AutoMPG',... % 1
       'bank',...% 2
       'diabetes',... % 3
       'Istanbul Stock Exchange',... % 4
       };
DataNameId      = 2;  
nn.DataNameId   = DataNameId;

% 样本数为偶数 
% 1= 'sinc', 2='AutoMPG', 3='bank', 4= 'diabetes',5='triazines' 
% 6='NonLinearSysIdentify', 7= 'Sincxy'   
DataName = string(DName(DataNameId))
Iter=1000;
RelAll=zeros(Iter,4);
for iter=1:Iter;
switch DataName
    case  'AutoMPG'  % 1  392x8
        trainNum=300; 
        AutoMPG=load('auto-mpg.data')';
        nn.mapmmPSFlag{DataNameId}=1;
        rand_sequence=randperm(size(AutoMPG,2),trainNum);
        rand_seqtest=setdiff(1:size(AutoMPG,2),rand_sequence);
        [traindata,nn.traindata_PS{DataNameId}] = mapminmax(AutoMPG(1:end-1,rand_sequence),0,1);
        [trainlabel,nn.trainlabel_PS{DataNameId}] = mapminmax(AutoMPG(end,rand_sequence),0,1);
        [testdata,nn.testdata_PS{DataNameId}] = mapminmax(AutoMPG(1:end-1,rand_seqtest),0,1);
        [testlabel,nn.testlabel_PS{DataNameId}] = mapminmax(AutoMPG(end,rand_seqtest),0,1); 
        J_max=17;
        Iter=5;
        % RSVD
        [Resultsrsvd,Net]=IT2TSK_NeuroFM_RSVD_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrersvd=PredictIT2_RSVD(size(testdata,2),testdata',testlabel',Net);
        J_max=17;
        Iter=5;
        [Resultsls,Net]=IT2TSK_NeuroFM_LS_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrels=PredictIT2_LS(size(testdata,2),testdata',testlabel',Net);
        iid=2; %=1 MSE,2: RMSE
        %                         
        RelAll(iter,:)=sqrt([Resultsrsvd(iid), ResultsPrersvd(iid),...% RSVD:Training (RMSE)   Testing (RMSE)
                             Resultsls(iid),ResultsPrels(iid)]);% LS:Training (RMSE)   Testing (RMSE)
    case  'bank' % 2
        bank=load('bank.data');  %8192x8
        bankIdent1=bank';
        trainNum=5E3;
        bankSize=size(bankIdent1,2);        
        rand_sequence=randperm(bankSize,trainNum);
        rand_seqtest=setdiff(1:bankSize,rand_sequence);
        nn.mapmmPSFlag{DataNameId}=1;
        [traindata,nn.traindata_PS{DataNameId}] = mapminmax(bankIdent1(1:end-1,rand_sequence),0,1);
        [trainlabel,nn.trainlabel_PS{DataNameId}] = mapminmax(bankIdent1(end,rand_sequence),0,1);
        [testdata,nn.testdata_PS{DataNameId}] = mapminmax(bankIdent1(1:end-1,rand_seqtest),0,1);
        [testlabel,nn.testlabel_PS{DataNameId}] = mapminmax(bankIdent1(end,rand_seqtest),0,1);  
        J_max=11;
        Iter=5;
        [Resultsrsvd,Net]=IT2TSK_NeuroFM_RSVD_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrersvd=PredictIT2_RSVD(size(testdata,2),testdata',testlabel',Net)
        J_max=11;
        Iter=5;
        [Resultsls,Net]=IT2TSK_NeuroFM_LS_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrels=PredictIT2_LS(size(testdata,2),testdata',testlabel',Net)
        iid=2; %=1 MSE,2: RMSE
        %                         
        RelAll(iter,:)=sqrt([Resultsrsvd(iid), ResultsPrersvd(iid),...% RSVD:Training (RMSE)   Testing (RMSE)
                             Resultsls(iid),ResultsPrels(iid)]);% LS:Training (RMSE)   Testing (RMSE)
    case 'diabetes' % 3  
        % diabetes
        diabetes2_data;     %   randomly generate new training and testing data for every trial of simulation
        traindata1=load('diabetes_train')'; %576x8
        testdata1=load('diabetes_test')';   %192x8
        nn.mapmmPSFlag{DataNameId}=0;
        traindata=traindata1(1:8,:);
        trainlabel=traindata1(end,:);
        testdata=testdata1(1:8,:);
        testlabel=testdata1(end,:); 
        J_max=11;
        Iter=5;
        [Resultsrsvd,Net]=IT2TSK_NeuroFM_RSVD_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrersvd=PredictIT2_RSVD(size(testdata,2),testdata',testlabel',Net);
        J_max=17;
        Iter=5;
        [Resultsls,Net]=IT2TSK_NeuroFM_LS_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrels=PredictIT2_LS(size(testdata,2),testdata',testlabel',Net)
         iid=2; 
        %%=1 MSE,2: RMSE                   
        RelAll(iter,:)=sqrt([Resultsrsvd(iid), ResultsPrersvd(iid),...% RSVD:Training (RMSE)   Testing (RMSE)
                             Resultsls(iid),ResultsPrels(iid)]);% LS:Training (RMSE)   Testing (RMSE)
    case 'Istanbul Stock Exchange' % 4:536x9
    % Attribute Information:
    % Stock exchange returns. Istanbul stock exchange national 100 index, Standard & poors 500 return index, 
    % Stock market return index of Germany, Stock market return index of UK, Stock market return index of Japan, 
    % Stock market return index of Brazil, MSCI European index, MSCI emerging markets index
        filename = 'ISTANBULSTOCK\data_akbilgic.xlsx';
        NonLinSysIdent = xlsread(filename)';
        trainNum=400;
        rand_sequence=randperm(size(NonLinSysIdent,2),trainNum);
        rand_seqtest=setdiff(1:size(NonLinSysIdent,2),rand_sequence);
        nn.mapmmPSFlag{DataNameId}=1;
        [traindata,nn.traindata_PS{DataNameId}] = mapminmax(NonLinSysIdent(2:end,rand_sequence),-1,1);
        [trainlabel,nn.trainlabel_PS{DataNameId}] = mapminmax(NonLinSysIdent(1,rand_sequence),-1,1);
        [testdata,nn.testdata_PS{DataNameId}] = mapminmax(NonLinSysIdent(2:end,rand_seqtest),-1,1);
        [testlabel,nn.testlabel_PS{DataNameId}] = mapminmax(NonLinSysIdent(1,rand_seqtest),-1,1);
        J_max=11;
        Iter=5;
        [Resultsrsvd,Net]=IT2TSK_NeuroFM_RSVD_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrersvd=PredictIT2_RSVD(size(testdata,2),testdata',testlabel',Net)
        J_max=11;
        Iter=5;
        [Resultsls,Net]=IT2TSK_NeuroFM_LS_testExample5func(J_max,Iter,traindata,trainlabel,nn);
        ResultsPrels=PredictIT2_LS(size(testdata,2),testdata',testlabel',Net)
        iid=2; %=1 MSE,2: RMSE
        %                         
        RelAll(iter,:)=sqrt([Resultsrsvd(iid), ResultsPrersvd(iid),...% RSVD:Training (RMSE)   Testing (RMSE)
                             Resultsls(iid),ResultsPrels(iid)]);% LS:Training (RMSE)   Testing (RMSE)
    otherwise
        warning('Wrong datasets, exit.')
        exit
end
end
format short E
min(RelAll)
%% write data to tex 
laxtabname2=strcat('Table5ID',num2str(DataNameId),'results.tex')
if exist(laxtabname2)
    delete(laxtabname2)
end
% Table structure:  L x noise_num
Results1m=mean(RelAll(:,1:2));%Example C RSVD
Results1s=std(RelAll(:,1:2));%Example C RSVD
Results2m=mean(RelAll(:,3:4));%Example C LS
Results2s=std(RelAll(:,3:4));%Example C LS   
str2  = strcat(    num2str(Results1m(1),'%.4d'),...
                  '   &', num2str(Results1s(1),'%.4d'),...
                  '   &', num2str(Results1m(2),'%.4d'),...
                  '   &', num2str(Results1s(2),'%.4d'),'\\');
dlmwrite(laxtabname2,str2,'delimiter', '','-append');
dlmwrite(laxtabname2, [],'newline','pc','-append');
str2  = strcat(    num2str(Results2m(1),'%.4d'),...
                      '   &', num2str(Results2s(1),'%.4d'),...
                      '   &', num2str(Results2m(2),'%.4d'),...
                      '   &', num2str(Results2s(2),'%.4d'),'\\');
    dlmwrite(laxtabname2,str2,'delimiter', '','-append');
dlmwrite(laxtabname2, [],'newline','pc','-append');
%dlmwrite(laxtabname2,'\bottomrule','delimiter', '','-append'); 
type(laxtabname2)


