sp=load('StockPrice.dat');
spsize=size(L);
Dir='StockPriceData';
if ~exist(Dir)
    mkdir(Dir);
end
tau=5;
for i=1:spsize(2) % ten  companies
    icomlen=length(sp(:,i));
    datai=[];
    for j=1:icomlen-tau
        datai(j,:)=[sp(j:j+tau,i)'];
    end
    filename=strcat(Dir,'\StockPriceData',num2str(i),'.mat');
    save(filename,'datai','sp_PS');
end
%% 
for i=1:spsize(2) % ten  companies
    filename=strcat(Dir,'\StockPriceData',num2str(i),'.mat');
    data=load(filename); 
end