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