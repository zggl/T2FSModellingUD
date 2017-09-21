function [X_scaled,Xij]=UniformDesignWithScale(n,s,coli,min_ranges_p,max_ranges_p);
uns = UDall(n,s);
Xij=(uns-0.5)/n;
Xij_i=Xij(:,:,coli);

%% scale to any interval 
% in_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the minimum value of each variable
% max_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the maximum value of each variable
% coli        : column index for uniform variable matrix.
% Xij         : X location cal from Uij 
% Latin Hypercube
% version 1.0 (2.38 KB) by Nassim Khaled
% lhsdesign_modified generates latin hypercube samples
% http://cn.mathworks.com/matlabcentral/fileexchange/45793-latin-hypercube
%  
plen=length(min_ranges_p);
if s~=plen
    error(spritf('minimum value of each variable should be %d',s) ) 
end
if s~=length(max_ranges_p)
    error(spritf('maximum value of each variable should be %d',s) ) 
end

[M,N]=size(min_ranges_p);
if M<N
    min_ranges_p=min_ranges_p';
end
    
[M,N]=size(max_ranges_p);
if M<N
    max_ranges_p=max_ranges_p';
end

slope=max_ranges_p-min_ranges_p;
offset=min_ranges_p;

SLOPE=ones(n,plen);
OFFSET=ones(n,plen);

for i=1:plen
    SLOPE(:,i)=ones(n,1).*slope(i);
    OFFSET(:,i)=ones(n,1).*offset(i);
end
% scaled results
X_scaled=SLOPE.*Xij_i+OFFSET;