%% predict
function [ResultsPre,Net]=PredictIT2_LS(N,In,Out,Net)
for ii=1:N
    x_t=In(ii,:)';
    if size(x_t,2)~=1
        error('[%s] must be column vector',num2str(x_t));
    end
    F.uk= prod( exp( -( (x_t*ones(1,Net(1).J_max)-Net(2).m{1}) ./ Net(2).sigma{1} ).^2 ) ,1); % lower
    F.lk= prod( exp( -( (x_t*ones(1,Net(1).J_max)-Net(2).m{2}) ./ Net(2).sigma{2} ).^2 ) ,1 ); % upper 
    fbar=Net(1).mm*F.uk/sum(F.uk)+Net(1).nn*F.lk/(sum(F.lk)+eps);
    A(ii,:)=kron(fbar,[1,x_t']);
end
%% Optimize $C_j$£¬ Recursive sigular value decomposition
Net(3).Oj3=A;
H0=A(1:Net(1).N0,:); 
T0=Out(1:Net(1).N0,:);
M = pinv(H0' * H0);
Net(3).LSYi{1}=H0;% for LS convergence ana
Net(3).LSXi{1}=T0;% for LS convergence ana
Net(3).Oj4=A*Net(3).Cj;
if Net(4).mapmmPSFlag{Net(4).DataNameId}==1
    Oj4=mapminmax('reverse',Net(3).Oj4,Net(4).trainlabel_PS{Net(4).DataNameId});
    Out=mapminmax('reverse',Out,Net(4).trainlabel_PS{Net(4).DataNameId});
    Net(4).ErrorSquareVec=(Oj4-Out).^2;
    Net(4).MSE=sum(Net(4).ErrorSquareVec)/(N);
    Net(4).RMSE=sqrt(sum(Net(4).ErrorSquareVec)/(N));
    Net(4).NRMSE=sqrt(sum(Net(4).ErrorSquareVec)/sum((Oj4-mean(Out)).^2));
    MSE=sum((Oj4-Out).^2)/(N);
else
    Net(4).ErrorSquareVec=(Net(3).Oj4-Out).^2;
    Net(4).MSE=sum(Net(4).Error)/(N);
    Net(4).RMSE=sqrt(sum(Net(4).Error)/(N));
    Net(4).NRMSE=sqrt(sum(Net(4).Error)/sum((Net(3).Oj4-mean(Out)).^2));
    MSE=Net(4).MSE;
end
Net(4).RMSE=sqrt(Net(4).MSE);
ResultsPre=[Net(4).MSE,Net(4).RMSE]; 
%disp('     MSE     RMSE')
end