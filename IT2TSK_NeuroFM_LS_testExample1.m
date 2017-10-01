%% SECTION TITLE
clc,clear
%% praparing data
J_max=8;
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
%% Type reduction method:
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
plot(X_scaled(:,1),X_scaled(:,2),'r*')
pause(1)
close all
%% 
Net(1).N0=55;    % Online Sequence LS
Net(1).Block=1;
Net(4).DataNameId=1;
Net(4).mapmmPSFlag{Net(4).DataNameId}=0;
Net=IT2TSKNeuroFM_LS(Net,QBFAStruc,X_scaled(:,1),Net(1).N0,Net(1).Block);
Net(2).m{1}
Net(2).sigma{1}
Net(3).Cj
sprintf('MSE: %4d, RMSE:  %4d, NRMSE: %4d ',Net(4).MSE,Net(4).RMSE,Net(4).NRMSE)
[ResultsPre,Net]=PredictIT2_LS(N,In',Out,Net);
ResultsPre
%% plot MFs

linsty={'-','--',':','-.','--',':','-.','-.',':','--'};
x=linspace(floor(QBFAStruc.bound(1)),ceil(QBFAStruc.bound(2)),100);

for inputnum=1:n
    figure(inputnum)
    legrul={};    
    for i=1:J_max
        if i<J_max
            if i==1
                legrul{i}=strcat(num2str(i),'st',' rule');
            elseif i==2
                legrul{i}=strcat(num2str(i),'nd',' rule');
            elseif i==3
                legrul{i}=strcat(num2str(i),'rd',' rule');
            else
                legrul{i}=strcat(num2str(i),'th',' rule');
            end
        else
            legrul{i}=strcat(num2str(J_max),'th',' rule');
        end
        y1=gaussmf(x,[Net(2).sigma{1}(inputnum,i) Net(2).m{1}(inputnum,i)]);
        y2=gaussmf(x,[Net(2).sigma{2}(inputnum,i) Net(2).m{1}(inputnum,i)]);
        plot([x,x(end:-1:1)],[y1,y2(end:-1:1)],'LineStyle',linsty{i})
        hold on
    end 
    legend(legrul)
    xlabel('$x$', 'Interpreter','latex' )
    ylabel('membership function value', 'Interpreter','latex' );
    epsname = strcat('UDSincMFs',num2str(inputnum),'.eps' );
    saveas(gcf,epsname,'epsc2') 
    hold off
end


