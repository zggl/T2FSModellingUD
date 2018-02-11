function [Pbest,Jbest]=QBFA_IT2TSKNN(QBFAStruc,Net)
clc;
EvalNum=QBFAStruc.EvalNum;
%IterMax                % The max iteration steps 
 for i=1:EvalNum        % eval it EvalNum(default is 50) times
     [Iter,Jbest,Pbest]=QBFA_main(QBFAStruc,Net,1);
end
end
%% Quantum inspired bacteria forage alogrithm
function [Iter,Jbest,Pbest]=QBFA_main(QBFAStruc,Net,id)
Iter=0;
C=0.01*rand(QBFAStruc.Nb,1);           % the run length
% bound of 
bound=QBFAStruc.bound;
DeflationSize=QBFAStruc.J_max*Net(1).DataDim(1);
BoundFun_m=(bound(2)-bound(1))*ones(1,DeflationSize);
BoundFun_sigma=(bound(4)-bound(3))*ones(1,DeflationSize);
BoundFun=[BoundFun_m,BoundFun_sigma];
BoundUB=[bound(2)*ones(1,DeflationSize),bound(4)*ones(1,DeflationSize)];
BoundLB=[bound(1)*ones(1,DeflationSize),bound(3)*ones(1,DeflationSize)];
% theta=linspace(0,pi/2,QBFAStruc.Nb);
J=zeros(QBFAStruc.Nb,QBFAStruc.Ncs+1);
P=zeros(QBFAStruc.Nb,QBFAStruc.VarNumber,QBFAStruc.Ncs+1);
if length(BoundFun)~=QBFAStruc.VarNumber
    error(sprintf('Center and width should be equal to %d',QBFAStruc.J_max));
end
% Population     the quantum individual
Population=zeros(2*QBFAStruc.Nb,QBFAStruc.VarNumber);
% Population, x,C,m,n,J
Jbest=Net(4).MSE;
Pbest=[Net(2).m{1}(:)',Net(2).sigma{1}(:)'];
for row=1:2:2*QBFAStruc.Nb
    if row==1
        Population(row,:)=Pbest;
        Population(row+1,:)=Pbest;
    else
        pop_i=Pbest+rand(1,QBFAStruc.VarNumber);
        Population(row,:)=pop_i;
        Population(row+1,:)=pop_i;
    end
end
% P(:,:,1)=[Pop1;
%           Population(row,:).^2.*(ones(length(row),1)*BoundFun)+bound(1)];
P(:,:,1)=Population(1:2:2*QBFAStruc.Nb,:);
Iter=Iter+QBFAStruc.Nb;
for eld=1:QBFAStruc.Ned
    for k=1:QBFAStruc.Nre
        for j=1:QBFAStruc.Ncs
            % U1----quantum rotate angle Eastern; U2----Clockwise
            theta1=0.005*pi;theta2=-0.005*rand*pi;
            U1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
            U2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
            if Iter>QBFAStruc.IterMax
               break;break;break;
            end
            for i=1:QBFAStruc.Nb
                %fprintf(1,'This is the  %d th bacteria''s %dth chemotaxis step  ( function id is:  %d ) \n',i,j,id);
%tumble
                J(i,j)=FN(P(i,:,j),QBFAStruc.J_max,Net);Iter=Iter+1;
                rowIndex=2*i-1:2*i;  
                if J(i,j)>Jbest
                    xLBestInd=find(P(i,:,j)<=Pbest);
                    xGBestInd=find(P(i,:,j)>Pbest);
                     if ~isempty(xLBestInd) %rotate Clockwise
                        Population(rowIndex,xLBestInd)=U2*Population(rowIndex,xLBestInd);
                     end
                     if ~isempty(xGBestInd) %rotate Eastern
                        Population(rowIndex,xGBestInd)=U1*Population(rowIndex,xGBestInd);
                     end
                else
                    xLBestInd=find(P(i,:,j)<=Pbest);
                    xGBestInd=find(P(i,:,j)>Pbest);
                     if ~isempty(xLBestInd)  %rotate Clockwise
                        Population(rowIndex,xLBestInd)=U1*Population(rowIndex,xLBestInd);
                     end
                     if ~isempty(xGBestInd) %rotate Eastern
                        Population(rowIndex,xGBestInd)=U2*Population(rowIndex,xGBestInd);
                     end                
                end
               P(i,:,j+1)=Population(rowIndex(1),:).^2.*BoundFun+bound(1);
               % check
               P(i,:,j+1)=PopulationChech(P(i,:,j+1),bound,BoundUB,BoundLB);
               [Pij,Net]=FN(P(i,:,j+1),QBFAStruc.J_max,Net);
               if Pij<Jbest
                   Jbest=Pij;
                   Pbest=P(i,:,j+1);   
%                    fprintf(1,'the  most best is (Jbest(tumble):%d), ( function id is:  %d ) \n',Jbest,id);
%                    fprintf(1,'the %dth bacteria,%dth chemotaxis step,%dth reproduction and %dth elimination step',i,j,k,eld);
               end
 %swim
          counter=0;Jlast=Pij;         % Initialize counter for swim length
            while counter<QBFAStruc.Ns
                  counter=counter+1;
                  if J(i,j+1)<=Jlast
                     Jlast=J(i,j+1);
                     Delta=(2*rand(1,QBFAStruc.VarNumber)-1).*rand(1,QBFAStruc.VarNumber);
                     P(i,:,j+1)=P(i,:,j+1)+C(i)*Delta/norm(Delta) ;
                     % check
                     P(i,:,j+1)=PopulationChech(P(i,:,j+1),bound,BoundUB,BoundLB);
                     [J(i,j+1),Net]=FN(P(i,:,j+1),QBFAStruc.J_max,Net);
                     Iter=Iter+1;
                     Jlast=J(i,j+1);Iter=Iter+1;
                     if Jlast<Jbest
                        Jbest=Jlast;
                        Pbest=P(i,:,j+1);                      
%                         fprintf(1,'the  most best is (Jbest(tumble):%d), ( function id is:  %d ) \n',Jbest,id);
%                         fprintf(1,'the %dth bacteria,%dth chemotaxis step,%dth reproduction and %dth elimination step',i,j,k,eld);
                     end
                  else
                     counter=QBFAStruc.Ns ;
                  end
            end
          end
%                 P(i,:,j)=Population(rowIndex(1),:).^2*BoundFun+bound(1);
%                 [Jbest,Jindex]=min(FN(P(:,:,j),id));Iter=Iter+QBFAStruc.Nb;
%                 Pbest=P(Jindex,:,j);
%                 P(i,:,j+1)=P(i,:,j);
        end
        P(:,:,1)=P(:,:,QBFAStruc.Ncs+1);
    end
%   Reprodution
        Jhealth=sum(J,2);              % Set the health of each of the S bacteria
        %C=C/1.1;
end
end
%%  Objective function
function [f_value,Net]=FN(varargin)
% x_t the ith point
% C   the cofficient of TSK rule
% Net  network structure
J=varargin{2};
Net=varargin{3};
% encoded parameters of IT2 MFs
MFPara_Nb=varargin{1};
f_Size=size(MFPara_Nb,1);
f_value=zeros(f_Size,1);
for jj=1:f_Size%
    MFPara=MFPara_Nb(jj,:);
    MFs_m=reshape(MFPara(1:Net(2).m_Len),Net(1).Dim{1},Net(1).Dim{3});
    MFs_sigma1=reshape(MFPara(1+Net(2).m_Len:end),Net(1).Dim{1},Net(1).Dim{3});
    MFs_sigma2=MFs_sigma1+Net(1).width;
    Net(2).MFs_m=MFs_m;
    Net(2).sigma1=MFs_sigma1;
    Net(2).sigma2=MFs_sigma2;
    for i=1:size(Net(1).Data{1},1)% for all the data  
        x_t=Net(1).Data{1}(i,:)';    
        %  BMM type reduction
        % Dongrui Wu; Mendel, J.M., "On the Continuity of Type-1 and Interval Type-2 Fuzzy Logic Systems," Fuzzy Systems, IEEE 
        if size(x_t,2)~=1
            error(sprintf('[%s] must be column vector',num2str(x_t)));
        end
        F.uk= prod( exp(-( (x_t*ones(1,J)-MFs_m) ./ MFs_sigma1 ).^2 ),1 ); % lower
        F.lk= prod( exp(-( (x_t*ones(1,J)-MFs_m) ./ (MFs_sigma2) ).^2 ),1 ); % upper 
        yk = Net(3).Cj*[1;x_t];
        y_t(i,1)=Net(1).mm*F.uk/sum(F.uk)*yk+Net(1).nn*F.lk/sum(F.lk)*yk;
    end % i
    f_value(jj)=sqrt(sum((Net(1).Data{2}-y_t).^2)/Net(1).Numb); 
end % jj
end
%% check population, single case
function x=PopulationChech(x,bound,BoundUB,BoundLB) 
    [pval_u,pinds_u] =find(x>BoundUB);
    [pval_l,pinds_l] =find(x<BoundLB);
    if ~isempty(pinds_u) 
        x(pinds_u) = 0.9* BoundUB(pinds_u); 
    end
    if ~isempty(pinds_l)
        x(pinds_l) =  1.2*BoundUB(pinds_l); 
    end    
end