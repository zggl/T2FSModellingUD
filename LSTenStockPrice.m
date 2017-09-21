clc
clear
TenStockP=load('StockPriceData/TenStockPriceErrorVec.mat','TenStockPri');
L=TenStockP.TenStockPri(1,:)';
[T,n]  = size(L);

N = 40; 
xls = zeros(T-N,n);
q = zeros(T-N,N);
qbar = zeros(T-N,N);
Q = zeros(T-N,1);
%Computing empirical frequencies
linsty={'-','--',':','-.','--',':','-.','-.',':','--'};
figure; 
hold on;
ten_ind=0;
step1=3;
for ten=1:step1:size(TenStockP.TenStockPri,1)
    ten_ind=ten_ind+1;
    legrul{ten_ind}=strcat(num2str(ten),'th',' stock prices');
    aux = zeros(2,N);
    L=TenStockP.TenStockPri(ten,:)';
    for i = 1:T-N
        LiN=L(i:i+N-1,:);
        xls(i,:) = sum(LiN,1)/N;
        q(i,:) = sort(diag((LiN-ones(N,1)*xls(i,:))*(LiN-ones(N,1)*xls(i,:))'));
        qbar(i,:) = ((N+5)/(N-1))*q(i,:); %See Example 1
        Q(i) = (L(i+N,:)-xls(i,:)) * (L(i+N,:)-xls(i,:))';
        aux = aux + [(Q(i)<=q(i,:)); (Q(i)<=qbar(i,:))]/(T-N);
    end 
    mcov = [aux(1,:); [1:1:N]/(N+1); aux(2,:)];
    xind=[1:N,N:-1:1];
    step2=4;
    Probx=8:step2:40;
    Proby(ten_ind,:)=mcov(3,Probx);
    % plot data in tables
    plot(xind,[mcov(1,:),mcov(3,end:-1:1)],'LineStyle',linsty{ten}) 
    if ten==1
        plot(mcov(2,:),'--b')
    end     
end
legend(legrul,'Location','northwest')
xlabel('$N$', 'Interpreter','latex' )
ylabel('Emprical costs', 'Interpreter','latex' );
box on
pause(3)
epsname1=strcat('StockPriceLSAsses','.eps' );
saveas(gcf,epsname1,'epsc2')
close all
%% write data to tex
probdata=[Probx;Proby;Probx/(N+1)];
laxtabname2='Table9results.tex'
if exist(laxtabname2)
    delete(laxtabname2)
end
dlmwrite(laxtabname2,probdata,'delimiter', '&');
 

