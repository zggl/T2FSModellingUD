%% SECTION TITLE
clc,clear
%% praparing data
J_max=4;
DNum=100;
In=linspace(0,2,DNum);
In_dim=size(In');
Out=zeros(In_dim);
n=In_dim(2);

for i=1:In_dim(1)
   etat=0.2*rand;
   Out(i,1)=0.9*In(i)/(1.2*In(i)^3+In(i)+0.3)+etat; 
end
N=DNum;
%surfc(X1,X2,Out1)
%% Type reduction method:
Struc.bound=[min(In(:));max(In(:));[0.3;0.8]];% bound of center m and width sigma
%% 
Net(1).IterTrainMax=10;
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
%% UD parameters
coli=1;
s=n+1;
if s==2
    Net(1).temp=minmax(In);
    min_ranges_p=Net(1).temp(:,1)'*ones(1,2);   
    max_ranges_p=Net(1).temp(:,2)'*ones(1,2);   
elseif s==3
    min_ranges_p=[-15,10,10]  
    max_ranges_p=[15, 30,20]  
end
[X_scaled,Xij]=UniformDesignWithScale(J_max,s,coli,min_ranges_p,max_ranges_p);
plot(X_scaled(:,1),X_scaled(:,2),'r*')
pause(1)
close all
%% 
N0=55;    % Online Sequence LS
Block=1;
Net(4).DataNameId=1;
Net(4).mapmmPSFlag{Net(4).DataNameId}=0;
Net=IT2TSKNeuroFM_LS(Net,Struc,X_scaled(:,1),N0,Block);

Net(2).m{1}
Net(2).sigma{1}
Net(3).Cj
sprintf('%4d',Net(4).MSE)
sprintf('%4d',Net(4).RMSE)
sprintf('%4d',Net(4).NRMSE)

%% plot MFs

linsty={'-','--',':','-.','--',':','-.','-.',':','--'};
x=linspace(floor(Struc.bound(1)),ceil(Struc.bound(2)),100);

for inputnum=1:n
    figure(inputnum)
    legrul={};    
    for i=1:J_max
        if i<J_max
            legrul{i}=strcat(num2str(i),'th',' rule');
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
    epsname = strcat('UDNonLinSysMFs',num2str(inputnum),'.eps' );
    saveas(gcf,epsname,'epsc2') 
    hold off
end
