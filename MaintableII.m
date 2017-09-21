clc,clear
laxtabname2='Table2results.tex'
if exist(laxtabname2)
    delete(laxtabname2)
end
    % Table structure:  L x noise_num
Jmax=[6:2:10];
for id=1:3 
    Results1=IT2TSK_NeuroFM_RSVD_testExample2func(Jmax(id),100);%Example 1 RSVD
    Results2=IT2TSK_NeuroFM_LS_testExample2func(Jmax(id),100);%Example 1 LS
    Results3=IT2TSK_NeuroFM_CWLS_testExample2func(Jmax(id),100);%Example 1 CWLS
    R1min1=min(Results1);
    R1min2=min(Results2);
    R1min3=min(Results3);    
    str2  = strcat(    num2str(R1min1(1),'%.4d'),...
                      '   &', num2str(R1min1(2),'%.4d'),...
                      '   &', num2str(R1min2(1),'%.4d'),...
                      '   &', num2str(R1min2(2),'%.4d'),...
                       '   &', num2str(R1min3(1),'%.4d'),...
                      '   &', num2str(R1min3(2),'%.4d'),'\\');
    dlmwrite(laxtabname2,str2,'delimiter', '','-append');
    dlmwrite(laxtabname2, [],'newline','pc','-append');
end
%dlmwrite(laxtabname2,'\bottomrule','delimiter', '','-append'); 