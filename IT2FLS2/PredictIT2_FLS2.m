%% predict
function [Net]=PredictIT2_RSVD(N,In,Out,Net)
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
% for t=1:N
%     Atemp=[];
%     inds_t=t;
%     if t==1
%         Intemp=A(inds_t,:);
%     else
%         Intemp=[Sugmati_t*V1';A(inds_t,:)];
%     end  
%    [U1,S1,V1]=svd((Intemp),'econ');
%    Sugmati_t=S1;  % \Sigma^{'}(t)
%    if t==1
%        Btemp=U1'*Out(t);
%    else
%        Btemp=U1'*[Btemp;Out(inds_t)];
%    end
% end
% C_star=V1*pinv(S1)*Btemp;
% Net(3).Cj=reshape(C_star',Net(1).J_max,[]);
if Net(4).mapmmPSFlag{Net(4).DataNameId}==1
    y_t=mapminmax('reverse',A*Net(3).Cj(:),Net(4).trainlabel_PS{Net(4).DataNameId});
    Out=mapminmax('reverse',Out,Net(4).trainlabel_PS{Net(4).DataNameId});
    Net(4).ErrorSquareVec=(y_t-Out).^2;
    Net(4).MSE=sum(Net(4).ErrorSquareVec)/N;
else
    y_t=A*Net(3).Cj(:);
    Net(4).ErrorSquareVec=(y_t-Out).^2;
    Net(4).MSE=sum(Net(4).ErrorSquareVec)/N;    
end
Net(4).RMSE=sqrt(Net(4).MSE);
%ResultsPre=[Net(4).MSE,Net(4).RMSE]; 
%disp('     MSE     RMSE')
end