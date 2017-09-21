J_max=7
Iter=2
Iteri=10;
Resplot=zeros(10,2);
ErrorRec=Resplot;
Net={};
for i=1:Iteri
    [Resultsi,datai_PS,N,Neti]=IT2TSK_NeuroFM_LS_StockPricefunc(J_max,Iter,i);
    Resplot(i,:)=min(Resultsi(:,:,i));
    ErrorRec(i,:)=[Neti(4).MSE,Neti(4).NRMSE];
    Net{i,1}=mapminmax('reverse',Neti(1).Data{2},datai_PS);
    Net{i,2}=mapminmax('reverse',Neti(3).Oj4,datai_PS);
end
% c = rand(Iteri,3);       %随机生成了12种颜色。RGB随机。
% c=[255,0,0;
%     255,127,80;%Coral 	珊瑚
%     255,140,0;% 	DarkOrange
%     218,165,32 % 	GoldEnrod秋麒麟
%     128,128,0; %  Olive橄榄 #808000 
%     127,255,0;%Chartreuse 	查特酒绿 	#7FFF00 
%     34,139,34;% ForestGreen 	森林绿 	#228B22 	
%     47,79,79;%DarkSlateGray 	深石板灰 	#2F4F4F 	
%     0,191,255;%DeepSkyBlue 	深天蓝 	#00BFFF 	
%     65,105,225;%RoyalBlue 	皇军蓝 	#4169E1 	
%     ]/255
dataori=[Net{:,1}];
datapre=[Net{:,2}];
plot(dataori,':'); hold on 
plot(datapre,'-'); hold off
xlim([0,N])
epsname1=strcat('StockPricePredicLS','.eps' );
saveas(gcf,epsname1,'epsc2')

%% write table
dlmwrite('StockPriceMSELS.tex',ErrorRec','delimiter','&')


