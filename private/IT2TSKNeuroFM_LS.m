function [Net]=IT2TSKNeuroFM_LS(varargin)
%     ocnzhao@gmail.com
%     12/06/2014
Net=varargin{1};
m_UD=varargin{3};
N0=varargin{4};
Block=varargin{5};
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
                        Net(1).Sj(Gt_ind)=St+1;
                        Net(2).sigma{1}(:,J)=Sigma0+sqrt( ( (St-1)*(sigma_t-Sigma0).^2+St*...
                            m_t2.^2+pv.^2 )/St-(St+1) / St * ( (St*m_t2+pv) / (St+1) ).^2 );
                        Net(2).sigma{2}(:,J)=Sigma0;
                    else % failes the output similarity test, greater or equal than tau*d
                        [Gjpv,Gj_ind]=min(abs(Out_i-Net(3).Cj(1:J)));
                        St=Net(1).Sj(Gj_ind); 
                        Net(1).Sj(Gj_ind)=St+1;
                    end        
                else % failes the input similarity test; less than rho
                    if J < J_max %  generate a new rule
                        J=J+1;
                        Net(2).sigma{1}(:,J)=Sigma0;
                        Net(2).sigma{2}(:,J)=Sigma0+Net(1).width;
                        Net(1).Sj(J)=1;
                    else % modify the Gt_ind th cluster 
                        St=Net(1).Sj(Gt_ind);
                        Net(2).sigma{2}(:,Gt_ind)=Sigma0;
                        Net(2).sigma{1}(:,Gt_ind)=Sigma0+Net(1).width;
                        Net(1).Sj(Gt_ind)=St+1; 
                    end  %%% J < J_max     
                end %%% ~isempty(indGj)      
            end  %%%% J==0   
        end %%%% ii
    end
%     Net(2).m{1}
%     Net(2).sigma{1}
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
%%  Optimize $m_ij,sigma_ij$
for IterTrain=1:Net(1).IterTrainMax
    Net(4).DataNameId=1;
    Net(4).mapmmPSFlag{Net(4).DataNameId}=0;
    inds=randperm(N);
    [MSE,Net] = LS(N,In,inds,Out,J,J_max,In_dim,N0,Block,Net);
end  % IterTrain
end 

%% RVD-Based Estimation
function [MSE,Net]=LS(N,In,inds,Out,J,J_max,In_dim,N0,Block,Net)
%% Optimize $C_j$£¬ Rrecursive SVD-Based Estimation
for ii=1:N
    x_t=In(inds(ii),:)';
    if size(x_t,2)~=1
        error('[%s] must be column vector',num2str(x_t));
    end
    F.uk= prod( exp( -( (x_t*ones(1,J_max)-Net(2).m{1}) ./ Net(2).sigma{1} ).^2 ) ,1); % lower
    F.lk= prod( exp( -( (x_t*ones(1,J_max)-Net(2).m{2}) ./ Net(2).sigma{2} ).^2 ) ,1 ); % upper 
    fbar=Net(1).mm*F.uk/sum(F.uk)+Net(1).nn*F.lk/sum(F.lk);
    A(inds(ii),:)=kron(fbar,[1,x_t']);
end
if rank(A)<J_max*(In_dim(2)+1)
    error('rank of A is inefficient, < %d',J_max*(In_dim(2)+1));
end
%% Optimize $C_j$£¬ Recursive sigular value decomposition
Net(3).Oj3=A;
H0=A(inds(1:N0),:); 
T0=Out(inds(1:N0),:);
M = pinv(H0' * H0);
beta = pinv(H0) * T0;
Net(3).LSXi{1}=H0;% for LS convergence ana
Net(3).LSYi{1}=T0;% for LS convergence ana
block_ind=1;
for n = N0 : Block : N
    if (n+Block-1) > N
        H = A(inds(n:N),:);    
        Tn = Out(inds(n:N),:);
        Block = size(H,1);             %%%% correct the block size 
    else
        H = A(inds(n:(n+Block-1)),:);    
        Tn = Out(inds(n:(n+Block-1)),:);
    end
    block_ind=1+block_ind;
    Net(3).LSXi{block_ind}=H;% for LS convergence ana
    Net(3).LSYi{block_ind}=Tn;% for LS convergence ana   
    M = M - M * H' * (eye(Block) + H * M * H')^(-1) * H * M; 
    beta = beta + M * H' * (Tn - H * beta);
end
% Kl=zeros(length(beta));
% for ik=1:block_ind
%     Kl=Kl+Net(3).LSXi{block_ind}'*Net(3).LSXi{block_ind}; 
%     qi(ik,1)=sum((Net(3).LSYi{block_ind}-Net(3).LSXi{block_ind}*beta).^2);
% end
% for ik=1:block_ind
%     vi=Net(3).LSXi{block_ind}'*Net(3).LSYi{block_ind};% dx1
%     hi=sum((Net(3).LSYi{block_ind}-Net(3).LSXi{block_ind}*vi).^2);%h_i:1x1
%     Ki=Net(3).LSXi{block_ind}'*Net(3).LSXi{block_ind};% dxd
%     Kibar=Ki+6*Ki*((Kl-Ki)\Ki);
%     if eig(Kl/6-Ki)>0
%         qibar(ik,1)=(beta-vi)'*Kibar*(beta-vi)+hi;
%     else                                        t
%         qiba      vbvbvbdfw r(ik,1)=inf;
%     end
%     Net(3).LSYi{block_ind}=H;% for LS convergence ana
%     Net(3).LSXi{block_ind}=Tn;% for LS convergence ana   
% end
Net(3).Oj4=A*beta;
if Net(4).mapmmPSFlag{Net(4).DataNameId}==1
    Oj4=mapminmax('reverse',Net(3).Oj4,Net(4).trainlabel_PS{Net(4).DataNameId});
    Out=mapminmax('reverse',Out,Net(4).trainlabel_PS{Net(4).DataNameId});
    Net(4).ErrorSquareVec=(Oj4-Out).^2;
    Net(4).MSE=sum(Net(4).ErrorSquareVec)/(N);
    Net(4).RMSE=sqrt(sum(Net(4).ErrorSquareVec)/(N));
    Net(4).NRMSE=sqrt(sum(Net(4).ErrorSquareVec)/sum((Oj4-mean(Out)).^2));
    MSE=sum(Net(4).ErrorSquareVec)/(N);
else
    Net(4).ErrorSquareVec=(Net(3).Oj4-Out).^2;
    Net(4).MSE=sum(Net(4).ErrorSquareVec)/(N);
    Net(4).RMSE=sqrt(sum(Net(4).ErrorSquareVec)/(N));
    Net(4).NRMSE=sqrt(sum(Net(4).ErrorSquareVec)/sum((Net(3).Oj4-mean(Out)).^2));
    MSE=Net(4).MSE;
end
Net(3).Cj=beta;
end

