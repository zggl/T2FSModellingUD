function uns = UDall(n,s)
% 均匀分布表 好格子点法
% n 水平数 s 因素数
h = 1:n;
ind = find(gcd(h,n)==1); %寻找比 n 小且与 n 互素的数 h
hm = h(ind); %生成向量
m = length(hm); %生成向量维数, 不小于因素数 s
udt = mod(h'*hm,n); %好格子点法
ind0 = find(udt==0);
udt(ind0)=n; %生成均匀设计表 U(n^s)
%udt(end,:)=[ ]; %生成均匀设计表 U^*((n- 1)^s)
if s>m
    disp('s 必须小于或等于 m');
return;
else
    mcs =nchoosek(m,s);
    if mcs<1e5
        tind = nchoosek(1:m,s);
        [p,q] = size(tind);cd2 = zeros(p,1);
        for k=1:p
            UT = udt(1:n,tind(k,:));
            cd2(k,1) = UDCD2(UT);
        end
        tc=tind(find(abs(cd2 - min(cd2))<1e-5),:);
        for r=1:size(tc,1);
            uns (:,:,r)= udt(:,tc(r,:));
        end
        else
        for k = 1:n
            a = k;
            UT = mod(h'*a.^(0:s- 1),n);
            cd2(k,1) = UDCD2(UT);
        end
        tc = find(abs(cd2 - min(cd2))<1e-5);
        for r=1:size(tc,1);
            uns (:,:,r)= mod(h'*tc(r).^(0:s- 1),n);
        end
        ind0 = find(uns==0); uns(ind0)=n;
    end
end