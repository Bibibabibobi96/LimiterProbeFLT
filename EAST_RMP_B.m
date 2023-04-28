% clear
% shot=52340;

% tpoint = 3.3;
% global IU IL

[RMP.t_Ip,RMP.Ip]=readfromEAST(shot,'pcrl01', 'pcs_east');
RMP.twin=[0 max(RMP.t_Ip)];


%% -----------get RMP current---------- 
disp('strat getting RMP current data')
for i=1:8
    RMP.sig=strcat('IRMPL',num2str(i));
    [RMP.t_RMP,x]=readfromEAST(shot,RMP.sig,'east_1', RMP.twin);
    RMP.RMPL(i,:)=smooth(x,10);
    
    sig=strcat('IRMPU',num2str(i));
    [RMP.t_RMP,x]=readfromEAST(shot,RMP.sig, 'east_1', RMP.twin);
    RMP.RMPU(i,:)=smooth(x,10);
end

%% -------------------------------------

[~,Index] = min(abs(RMP.t_RMP-tpoint));
% global IU IL
RMP.IU = RMP.RMPU(:,Index);
RMP.IL = RMP.RMPL(:,Index);


subplot(2,1,1)
plot(RMP.t_RMP,RMP.RMPL)
title('Lower')
subplot(2,1,2)
plot(RMP.t_RMP,RMP.RMPL)
title('Uper')

%% -----------------------------------

% EAST_RMP_B

load('G:\LimiterProbeFLT\B_field\RMP_basic.mat')
% load('G:\LimiterProbeFLT\B_field\matlab.mat')
% load('G:\LimiterProbeFLT\B_field\RMP_zero.mat')

RMP.II = [RMP.IU;RMP.IL];

BR3D = zeros(size(BR3D_RMP{i},1),size(BR3D_RMP{i},2),size(BR3D_RMP{i},3));
BZ3D = BR3D;
Bphi3D = BR3D;

for i = 1:16
    BR3D = BR3D + RMP.II(i) * BR3D_RMP{i};
    BZ3D = BZ3D + RMP.II(i) * BZ3D_RMP{i};
    Bphi3D = Bphi3D + RMP.II(i) * Bphi3D_RMP{i};
end

% clear BR3D_RMP BZ3D_RMP Bphi3D_RMP sig

% save (['G:\ÏÞÖÆÆ÷\LimiterProbeCode\LimiterProbeFLT\Bfield\',filename2],'Bphi3D','BR3D','BZ3D')