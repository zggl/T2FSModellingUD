clc
clear

try
    cd SinCxy
end

li=pi;

lival=-li;rival=li;
%% training data
dataNumdim=41;
dataNumdimxy=dataNumdim*dataNumdim;
x=linspace(lival,rival,dataNumdim)';
y=x;
[X,Y] = meshgrid(x,y);
sincx=@(x)sin(x)/x;
Sincxytr=zeros(dataNumdimxy,3);
for i=1:dataNumdim
    for j=1:dataNumdim
        k=dataNumdim*(i-1)+j;
        Sincxytr(k,:)=[x(i),y(j),abs(sinc(pi*x(i))*sinc(pi*y(j)))];
    end
end
%surf(X,Y,reshape(Sincxytr(:,3),dataNumdim,[]),'FaceAlpha',0.75,'EdgeColor','flat')
%% testing data
lival=-li;rival=li;
dataNumdim=30;
dataNumdimxy=dataNumdim*dataNumdim;
x=linspace(lival,rival,dataNumdim)';
y=x;
Sincxyte=zeros(dataNumdimxy,3);
for i=1:dataNumdim
    for j=1:dataNumdim-1
        k=dataNumdim*(i-1)+j;
        Sincxyte(k,:)=[x(i),y(j),abs(sinc(pi*x(i))*sinc(pi*y(j)))];
    end
end
%surf(X,Y,reshape(Sincxyte(:,3),dataNumdim,[]),'FaceAlpha',0.75,'EdgeColor','flat')
save('SincxyData','Sincxytr','Sincxyte')
try
    cd ..
end

