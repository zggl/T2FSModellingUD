function [Net]=IT2TSKNeuroFM_RSVD(varargin)
%     ocnzhao@gmail.com
%     12/06/2014
Net=varargin{1};
m_UD=varargin{3};
J_max=Net(1).J_max;
In=Net(1).Data{1};
Out1=Net(1).Data{2};
%surfc(X1,X2,Out1)
In_dim=Net(1).In_dim;
rho=Net(1).rho;           % threshold
tau=Net(1).tau;           % threshold
N=Net(1).Numb;
Out=Out1(:);
d=max(Out)-min(Out);
Sigma0=Net(2).sigma0;
% Generate rough fuzzy rules
% synamic size of cluster number
J=0;  
redoflag=0;
Net(2).m{1}=m_UD';
Net(2).m{2}=m_UD';
while redoflag==0
    for IterLearning=1:5 
        inds=randperm(N);
        for ii=1:N
            inds_ii=inds(ii);
            % The i th training data
            In_i=In(inds_ii,:)';
            Out_i=Out(inds_ii)';
            if J==0 % The first  round
                J=J+1;
                Net(2).sigma{2}(:,J)=Sigma0;
                Net(2).sigma{1}(:,J)=Sigma0+Net(1).width;
                Net(3).Cj(J,:)=[1,0.1*ones(1,In_dim(2))]; 
                Net(1).Sj(J)=1;
            else                % clusters number >=1
                % Input similarity test
                GjPv=prod( exp(-( (In_i*ones(1,J)-Net(2).m{2}(:,1:J)) ./ Net(2).sigma{2}(:,1:J) ).^2 ) );
                Net(1).Gj(inds_ii,1:J)=GjPv;        % store the iith firing strength
                indGj=find(GjPv>=rho);
                [Gtpv,Gt_ind]=max(GjPv);            % modify the cluster Gt_ind
                if ~isempty(indGj)% passes the input similarity test
                    % Output similarity test
                    e_vj=abs(Out_i-Net(3).Cj(1:J));
                    [row,col] = find(e_vj<=tau*d);
                    if ~isempty(col)                % passes the output similarity test
                        St=Net(1).Sj(Gt_ind);
                        sigma_t=Net(2).sigma{2}(:,Gt_ind);
                        m_t1=Net(2).m{1}(:,Gt_ind);
                        m_t2=Net(2).m{2}(:,Gt_ind);
                        pv=In_i;
                        Net(3).Cj(Gt_ind,:)=[1,(Net(3).Cj(Gt_ind,2:end)*St+0.01*ones(1,In_dim(2)))/(St+1)];
                        Net(1).Sj(Gt_ind)=St+1;
                        Net(2).sigma{1}(:,J)=Sigma0+sqrt( ( (St-1)*(sigma_t-Sigma0).^2+St*...
                            m_t2.^2+pv.^2 )/St-(St+1) / St * ( (St*m_t2+pv) / (St+1) ).^2 );
                        Net(2).sigma{2}(:,J)=Sigma0;
                    else % failes the output similarity test, greater or equal than tau*d
                        [Gjpv,Gj_ind]=min(abs(Out_i-Net(3).Cj(1:J)));
                        St=Net(1).Sj(Gj_ind);
                        Net(3).Cj(Gj_ind,:)=[1,(Net(3).Cj(Gj_ind,2:end)*St+0.01*ones(1,In_dim(2)))/(St+1)];  
                        Net(1).Sj(Gj_ind)=St+1;
                    end        
                else % failes the input similarity test; less than rho
                    if J < J_max %  generate a new rule
                        J=J+1;
                        Net(2).sigma{1}(:,J)=Sigma0;
                        Net(2).sigma{2}(:,J)=Sigma0+Net(1).width;
                        Net(3).Cj(J,:)=[1,0.1*Out_i*ones(1,In_dim(2))]; 
                        Net(1).Sj(J)=1;
                    else % modify the Gt_ind th cluster 
                        St=Net(1).Sj(Gt_ind);
                        Net(2).sigma{2}(:,Gt_ind)=Sigma0;
                        Net(2).sigma{1}(:,Gt_ind)=Sigma0+Net(1).width;
                        Net(3).Cj(Gt_ind,:)=[1,(Net(3).Cj(Gt_ind,2:end)*St+0.01*ones(1,In_dim(2)))/(St+1)];
                        Net(1).Sj(Gt_ind)=St+1; 
                    end  %%% J < J_max     
                end %%% ~isempty(indGj)      
            end  %%%% J==0   
        end %%%% ii
    end
    if J<J_max
        %error('There should be %d Rules, rerun the script',J_max);
        redoflag=0;
    else
       redoflag=1; 
    end
    if IterLearning==5 
        break
    end
end
% Net(2).m{1}
% Net(2).sigma{1}
%%  Optimize $m_ij,sigma_ij$, QBFA
for IterTrain=1:Net(1).IterTrainMax
    [MSE,Net] = RSVD(N,In,Out,J,J_max,In_dim,Net);
    %MSEVec(IterTrain)=MSE;
end  % IterTrain
end 
% Net(2).m,Net(2).sigma,Net(3).Cj,Net(1).Sj,sum(Net(1).Sj);

%% RVD-Based Estimation
function [MSE,Net]=RSVD(N,In,Out,J,J_max,In_dim,Net)
%% Optimize $C_j$£¬ Rrecursive SVD-Based Estimation
for ii=1:N
    x_t=In(ii,:)';
    if size(x_t,2)~=1
        error('[%s] must be column vector',num2str(x_t));
    end
    F.uk= prod( exp( -( (x_t*ones(1,J_max)-Net(2).m{1}) ./ Net(2).sigma{1} ).^2 ) ,1); % lower
    F.lk= prod( exp( -( (x_t*ones(1,J_max)-Net(2).m{2}) ./ Net(2).sigma{2} ).^2 ) ,1 ); % upper 
    fbar=Net(1).mm*F.uk/sum(F.uk)+Net(1).nn*F.lk/(sum(F.lk)+eps);
    A(ii,:)=kron(fbar,[1,x_t']);
end
if rank(A)<J_max*(In_dim(2)+1)
    error('rank of A is inefficient, < %d',J_max*(In_dim(2)+1));
end
%% Optimize $C_j$£¬ Recursive sigular value decomposition
Net(3).Oj3=A;
for t=1:N
    Atemp=[];
    inds_t=t;
    if t==1
        Intemp=A(inds_t,:);
    else
        Intemp=[Sugmati_t*V1';A(inds_t,:)];
    end  
   [U1,S1,V1]=svd((Intemp),'econ');
   Sugmati_t=S1;  % \Sigma^{'}(t)
   if t==1
       Btemp=U1'*Out(t);
   else
       Btemp=U1'*[Btemp;Out(inds_t)];
   end
end
C_star=V1*pinv(S1)*Btemp;
Net(3).Cj=reshape(C_star',J_max,[]);
Net(3).Oj4=A*C_star;
if Net(4).mapmmPSFlag{Net(4).DataNameId}==1
    Oj4=mapminmax('reverse',Net(3).Oj4,Net(4).trainlabel_PS{Net(4).DataNameId});
    Out=mapminmax('reverse',Out,Net(4).trainlabel_PS{Net(4).DataNameId});
    Net(4).MSE=sum((Oj4-Out).^2)/(N);
    Net(4).RMSE=sqrt(sum((Oj4-Out).^2)/(N));
    Net(4).NRMSE=sqrt(sum((Oj4-Out).^2)/sum((Net(3).Oj4-mean(Out)).^2));
    Net(4).Error=Oj4-Out;
    MSE=sum(Net(4).Error.^2)/(N);
else
    Net(4).MSE=sum((Net(3).Oj4-Out).^2)/(N);
    Net(4).RMSE=sqrt(sum((Net(3).Oj4-Out).^2)/(N));
    Net(4).NRMSE=sqrt(sum((Net(3).Oj4-Out).^2)/sum((Net(3).Oj4-mean(Out)).^2));
    Net(4).Error=Net(3).Oj4-Out;
    MSE=sum(Net(4).Error.^2)/(N);   
end
end

