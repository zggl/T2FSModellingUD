clc
clear
DName={'AutoMPG',...    % 1
       'bank',...       % 2
       'diabetes',...   % 3
       'Istanbul Stock Exchange',... % 4
       };
DataNameId  = 2;        % change DataNameId to the datasets you want
name=strcat('Errorvector',num2str(DataNameId),'WLS');
epsname=strcat(name,'.mat');
Coverage=load(epsname,'AA'); % Nx1E4
CoverageAA=Coverage.AA;
if DataNameId == 1
    N = 39;
elseif DataNameId == 2
    N = 3099;
elseif DataNameId == 3
    N = 139;
elseif DataNameId == 4
    N = 99;    
end
if DataNameId == 1
    BINS=20;
    HEIGHT=15;
elseif DataNameId == 2
    BINS=30;
    HEIGHT=15;
elseif DataNameId == 3
    BINS=20;
    HEIGHT=10;   
elseif DataNameId == 4
    BINS=15; 
    HEIGHT=10;
end
if DataNameId == 1
    covprob=0.8;
elseif DataNameId == 2
    covprob=0.8;
elseif DataNameId == 3
    covprob=0.8;
elseif DataNameId == 4
    covprob=0.8;
end
Probx=floor((N+1)*covprob);
Iteraize=size(CoverageAA,2);
Proby=zeros(Iteraize,2);
for k=1:Iteraize
    L=CoverageAA(:,k);
    [T,n]  = size(L);
    xls = zeros(T-N,n);
    q = zeros(T-N,N);
    qbar = zeros(T-N,N);
    Q = zeros(T-N,1);
    %Computing empirical frequencies
    aux = zeros(2,N);
    for i = 1:T-N
        LiN=L(i:i+N-1,:);
        xls(i,:) = sum(LiN,1)/N;
        q(i,:) = sort(diag((LiN-ones(N,1)*xls(i,:))*(LiN-ones(N,1)*xls(i,:))'));
        qbar(i,:) = ((N+5)/(N-1))*q(i,:); %See Example 1
        Q(i) = (L(i+N,:)-xls(i,:)) * (L(i+N,:)-xls(i,:))';
        aux = aux + [(Q(i)<=q(i,:)); (Q(i)<=qbar(i,:))]/(T-N);
    end 
    mcov = [aux(1,:); [1:1:N]/(N+1); aux(2,:)];
    Proby(k,1:2)=mcov([1,3],Probx);
end
maxheight=0;
q_flag=1;
[fx,x]=hist(Proby(:,1),linspace(min(Proby(:,q_flag)'),max(Proby(:,q_flag)'),BINS)); % fx:counts; x:centers
areafx=trapz(x,fx);
maxheight=max(maxheight,max(fx/areafx));
histog{1}=fx/areafx;
histog{2}=x;
h=bar(histog{2},histog{1},'barwidth',1);  
bar_setFaceAlpha(h,0.1);
set(h,'FaceColor', [1 0 0]);
hold on
q_flag=2;
[fx,x]=hist(Proby(:,2),linspace(min(Proby(:,q_flag)'),max(Proby(:,q_flag)'),BINS)); 
areafx=trapz(x,fx);
maxheight=max(maxheight,max(fx/areafx));
histog{1}=fx/areafx;
histog{2}=x;
h=bar(histog{2},histog{1},'barwidth',1);
% title(['$q_{',num2str(Probx),'}$','with $N$=',num2str(N)],'Interpreter','LaTex') 
% title(['N=',num2str(N),' $\bar{q}_{',num2str(Probx),'}$'],'Interpreter','LaTex')
bar_setFaceAlpha(h,0.7);
set(h,'FaceColor', [0.38,0.396,0.506]);
xlabel('Coverage') 
ylabel('Empirical frequency')
axis([0,1,0,HEIGHT])
epsname1=strcat(name,'.eps' );
saveas(gcf,epsname1,'epsc2')
pause(2)
close all
