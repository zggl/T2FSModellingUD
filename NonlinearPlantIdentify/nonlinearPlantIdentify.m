ykp1=0;
ykp2=0;
type=2;
dNum=200
y=zeros(1,dNum);
dataNPlant=zeros(dNum,4);
 
% type=2 
for i=1:dNum
    if type==1
        ukp1=-2 + 4*rand;
    else
        ukp1=sin(2*pi*i/25);
    end  
    y(i)=(ykp1*ykp2*(ykp1+2.5))/(1+ykp1^2+ykp2^2)+ukp1;
    dataNPlant(i,:)=[ykp1,ykp2,ukp1,y(i)];
    ykp1=y(i);
    ykp2=ykp1;
    
end
plot(y)
xlabel('$k$', 'Interpreter','latex' )
%ylabel('$\frac{\sin(x)}{x}$', 'Interpreter','latex' )
ylabel('$y_p(k)$', 'Interpreter','latex' );
epsname = strcat('nonlinearfun',num2str(type), '.eps' );
saveas(gcf,epsname,'epsc2')


save(strcat('NonlinearDataIdentify',num2str(type)),'dataNPlant')

