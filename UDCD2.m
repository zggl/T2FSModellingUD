function CD2 = UDCD2(UT)
%ÖÐÐÄ»¯ L2 Æ«²î CD2
U = UT;
[n,s] = size(U);
X = (2*U- 1)/(2*n);
cs1 = zeros(n,1);
for k=1:n
    CSP=1;
    for i=1:s 
        CSP = CSP*(2+abs(X(k,i)- 1/2)- (X(k,i)- 1/2)^2);
    end
    cs1(k) = CSP;
end
CS1=sum(cs1);
cs2 = zeros(n,n);
for k=1:n
    for l=1:n
        CSP=1;
        for i=1:s
        CSP = CSP*(1+1/2*abs(X(k,i)- 1/2)+1/2*abs(X(l,i)- 1/2)-...
                1/2*abs(X(k,i)- X(l,i)));
        end
        cs2(k,l) = CSP;
    end
end
CS2 = sum(sum(cs2));
CS = (13/12)^s - 2^(1- s)/n*CS1 + 1/(n^2)*CS2;
CD2 = sqrt(CS);