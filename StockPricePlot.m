sp=load('StockPrice.dat');
plot(sp)
xlim([0,950])
xlabel('Day $i$', 'Interpreter','latex' )
ylabel('Stock Price', 'Interpreter','latex' );
epsname = strcat('StockPriceData','.eps' );
saveas(gcf,epsname,'epsc2') 