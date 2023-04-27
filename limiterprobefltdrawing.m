clear
clc

tic
% B_total;
shot=103950;
tpoint=5.4;


% shot=shotnum;
% -----------
global nphi
nphi=720;

filename1=[num2str(shot),'_',num2str(tpoint*1000),'.0.mat'];
% disp(filename1)
filename2=[num2str(shot),'at',num2str(tpoint),'_RMP.mat'];
% disp(filename2)

% get EFIT magnetic field
if exist(['G:\ERGOS_HCFs\dig\notebooks\',filename1])~=0
    load(['G:\ERGOS_HCFs\dig\notebooks\',filename1])
    cprintf('text','efit field file ');
    cprintf([1,0.5,0],filename1);
    cprintf('text',' is loaded\n');
else 
    error('EFIT file does not exist')
end

% get RMP field
if exist(['G:\限制器\LimiterProbeCode\LimiterProbeFLT\Bfield\',filename2])~=0
    load(['G:\限制器\LimiterProbeCode\LimiterProbeFLT\Bfield\',filename2])
    cprintf('text','RMP field file ');
    cprintf([1,0.5,0],filename2);
    cprintf('text',' is loaded\n');
else 
    disp('RMP field file does not exist, start calculating')
    EAST_RMP_current;
    disp('RMP field calculation finished')
%     load(['G:\限制器\LimiterProbeCode\LimiterProbeFLT\Bfield\',filename2])
end

% check if the size of RMP filed match with EFIT field
if length(BR3D(1,1,:)) ~= nphi
    error('nphi of 3D magnetic field and FLT msut be identical')
end

% modification of matrix in order of  R,Z,Phi
for i=1:nphi
    Bfield.BR3D_efit(:,:,i) = permute(BR, [2 1]);
    Bfield.BZ3D_efit(:,:,i) = permute(BZ, [2 1]);   
    Bfield.Bphi3D_efit(:,:,i) = permute(Bt, [2 1]);
end

% 4 turns of RMP coil
% modification of matrix in order of  R,Z,Phi
Bfield.BR3D_RMP = 4*permute(BR3D, [2 1 3]);
Bfield.BZ3D_RMP = 4*permute(BZ3D, [2 1 3]);
Bfield.Bphi3D_RMP = 4*permute(Bphi3D, [2 1 3]);

% integrate RMP B field and EFIT B field

Bfield.BR3D_total = Bfield.BR3D_RMP + Bfield.BR3D_efit;

Bfield.BZ3D_total = Bfield.BZ3D_RMP + Bfield.BZ3D_efit;

Bfield.Bphi3D_total = Bfield.Bphi3D_RMP + Bfield.Bphi3D_efit;


clear BR BZ Bt BR3D BZ3D Bphi3D 


%% -------------
% sectional view of B field 

figure(222)

i = 5
aa = nphi/2;
t = tiledlayout(3,3);

nexttile(1)
contourf(R,Z,Bfield.BR3D_efit(:,:,aa)','ShowText','on');
xlabel('R');ylabel('Z');title('B_{r}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.BZ3D_efit(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{z}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.Bphi3D_efit(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{t}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.BR3D_RMP(:,:,aa)','ShowText','on');
% contourf(R,Z,BR3D_RMP{i}(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{r}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.BZ3D_RMP(:,:,aa)','ShowText','on');
% contourf(R,Z,BZ3D_RMP{i}(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{z}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.Bphi3D_RMP(:,:,aa)','ShowText','on');
% contourf(R,Z,Bphi3D_RMP{i}(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{t}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.BR3D_total(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{r}');
colorbar
axis equal tight

nexttile
contourf(R,Z,Bfield.BZ3D_total(:,:,aa)');
xlabel('R');ylabel('Z');title('B_{z}');
colorbar
axis equal tight

nexttile
pcolor(R,Z,Bfield.Bphi3D_total(:,:,aa)'); shading interp
xlabel('R');ylabel('Z');title('B_{t}');
colorbar
axis equal tight

t.TileSpacing = 'compact';
t.Padding = 'none';


%% -----Interp3-----
% Interpolate field data to higher resolution

R_interp = linspace(min(R),max(R),129*5);
Z_interp = linspace(min(Z),max(Z),129*5);
[interp_r,interp_z,interp_t] =...
    meshgrid(linspace(min(R),max(R),129*5),linspace(min(Z),max(Z),129*5),linspace(0,2*pi,nphi));

BR3D_total = interp3(R,Z,linspace(0,2*pi,nphi),Bfield.BR3D_total,interp_r,interp_z,interp_t,'cubic');
BZ3D_total = interp3(R,Z,linspace(0,2*pi,nphi),Bfield.BZ3D_total,interp_r,interp_z,interp_t,'cubic');
Bphi3D_total = interp3(R,Z,linspace(0,2*pi,nphi),Bfield.Bphi3D_total,interp_r,interp_z,interp_t,'cubic');

clear interp_r interp_z interp_t
clear Bfield

BR3D_total(:,:,nphi+1) = BR3D_total(:,:,1);
BZ3D_total(:,:,nphi+1) = BZ3D_total(:,:,1);
Bphi3D_total(:,:,nphi+1) = Bphi3D_total(:,:,1);

%% Field Line Tracing
clear flt
tic
flt.n = 10;
flt.R_min = min(R);
flt.R_max = max(R);
flt.Z_min = min(Z);
flt.Z_max = max(Z);

% Define 3D bondary condition

% load('EAST_new_divertor.mat')
% load('EAST_new_divertor_limiter+0.01.mat')
% load('EAST_new_divertor_limiter+0.05_3D.mat')
% load('EAST_new_divertor_limiter+LHW+0.05_3D.mat')
% load('G:\LimiterProbeFLT\EAST_bondary\EAST_new_divertor_limiter+0.01_3D.mat')
load('G:\LimiterProbeFLT\EAST_bondary\EAST_square.mat')

flt.bondary_r = bondary.r;
flt.bondary_z = bondary.z;

% Define start points

% load ('limiterprobeposition.mat')
load('gogogo137.mat')
% load ('FLTtest.mat')
% load ('limiterprobeposition_inner0.1.mat')
flt.nr_mesh = 10;
flt.nz_mesh = 40;

flt.rStart = linspace(2.1, 2.35, flt.nr_mesh);
flt.zStart = linspace(-0.5, 0.5, flt.nz_mesh);
% flt.rStart = linspace(1.6, 1.7, flt.nr_mesh);
% flt.zStart = linspace(0, 0, flt.nz_mesh);


flt.tStart=-t+2*pi*(13/32);
%correction to real EAST toroidal angel

% t=t+2*pi*(13/32);
[flt.rStart,flt.zStart] = meshgrid(flt.rStart,flt.zStart);
flt.rStart_real = reshape(flt.rStart,1, size(flt.rStart,1) * size(flt.rStart,2));
flt.zStart_real = reshape(flt.zStart,1, size(flt.zStart,1) * size(flt.zStart,2));
flt.tStart_real = 0.0340 * ones(1,flt.nr_mesh * flt.nz_mesh);
flt.nphiStart = ceil(flt.tStart(1)/(2*pi/nphi));


clear h
global h
h=2*pi/nphi;

disp('Start Backward Tracing')
for i = 1 : flt.nr_mesh * flt.nz_mesh
    disp(['Start tracing from point',num2str(i)]);
    flt.rStart = flt.rStart_real(1,i);
    flt.zStart = flt.zStart_real(1,i);
    flt.tStart = flt.tStart_real(1,i);
    [flt.backward_r(:,i),flt.backward_z(:,i),flt.backward_t(:,i),flt.Lc(1,i)] = LPflt(flt,BR3D_total,BZ3D_total,Bphi3D_total);
%     toc
end
flt.backward_x = flt.backward_r .* cos(flt.backward_t);
flt.backward_y = flt.backward_r .* sin(flt.backward_t);

flt.Lc(2,:) = flt.Lc(1,:);

disp('Start Forward Tracing')
h=-h;
for i = 1 : flt.nr_mesh * flt.nz_mesh
    disp(['Start tracing from point',num2str(i)]);
    flt.rStart = flt.rStart_real(1,i);
    flt.zStart = flt.zStart_real(1,i);
    flt.tStart = flt.tStart_real(1,i);
    [flt.forward_r(:,i),flt.forward_z(:,i),flt.forward_t(:,i),flt.Lc(2,:)] = LPflt(flt,BR3D_total,BZ3D_total,Bphi3D_total);
end

flt.forward_x = flt.forward_r .* cos(flt.forward_t);
flt.forward_y = flt.forward_r .* sin(flt.forward_t);

flt.forward_x = flt.forward_x(end:-1:1,:);
flt.forward_y = flt.forward_y(end:-1:1,:);
flt.forward_z = flt.forward_z(end:-1:1,:);
flt.forward_r = flt.forward_r(end:-1:1,:);
flt.forward_t = flt.forward_t(end:-1:1,:);

flt.x = [flt.forward_x;flt.backward_x];
flt.y = [flt.forward_y;flt.backward_y];
flt.z = [flt.forward_z;flt.backward_z];
flt.r = [flt.forward_r;flt.backward_r];
flt.t = [flt.forward_t;flt.backward_t];


flt.Lc(3,:) = flt.Lc(1,:) + flt.Lc(2,:);
toc
%%

figure
plot(bondary.r(:,1),bondary.z(:,1));hold on
scatter(flt.backward_r(1:nphi:nphi*flt.n,:),flt.backward_z(1:nphi:nphi*flt.n,:),0.5);
scatter(flt.forward_r(1:nphi:nphi*flt.n,:),flt.forward_z(1:nphi:nphi*flt.n,:),0.5);

toc

%% ----------------magnetic field line plot-------------
disp(['Start magnetic field line plot']);
figure(111)

% for j=1:i
%     plot3(poincare_x(:,j),poincare_y(:,j),poincare_z(:,j))
% hold on
% end
plot3(flt.forward_x,flt.forward_y,flt.forward_z,'r');hold on
plot3(flt.backward_x,flt.backward_y,flt.backward_z,'b');
% plot3(flt_x,flt_y,flt_z,'r');
hold on

color = jet(n_point);
% color = hsv(n_point);

[limiterprobex,limiterprobey,limiterprobez] = pol2cart(t(1,:),r(1,:),z(1,:));
%plot limiter probe head
for i=1:length(r(1,:))
    plot3(limiterprobex(i),limiterprobey(i),limiterprobez(i),'o','MarkerFaceColor',color(i,:),'MarkerSize',8,'LineWidth',1.5,'MarkerEdgeColor','k');
    % scatter3(limiterprobex,limiterprobey,limiterprobez,s,c,'filled','MarkerEdgeColor','k')
    hold on
end

%%------------plot the end of field line---------
l=length(flt.x);
% Lcp;
% Lcn;

for i=1:length(r(1,:))
    plot3(flt.x(l/2-Lc(1,i)+1,i),flt.y(l/2-Lc(1,i)+1,i),flt.z(l/2-Lc(1,i)+1,i),'o','MarkerFaceColor',color(i,:),'MarkerSize',8,'LineWidth',1.5,'MarkerEdgeColor','k');% in positive direction
    hold on
end

for i=1:length(r(1,:))
    plot3(flt.x(l/2+Lc(2,i),i),flt.y(l/2+Lc(2,i),i),flt.z(l/2+Lc(2,i),i),'o','MarkerFaceColor',color(i,:),'MarkerSize',8,'LineWidth',1.5,'MarkerEdgeColor','k');% in positive direction
    hold on
end


%---------------------plot RMP,limiter and LHW antenna-------
load('RMPxyz.mat')
for j=1:16
    xyz=RMPxyz{j};
    plot3(xyz(1,:),xyz(2,:),xyz(3,:),'g','linewidth',4)
    hold on
end
for i = 1:size(bondary.LHW1_x,2)
    plot3(bondary.LHW1_x(:,i),bondary.LHW1_y(:,i),bondary.LHW1_z(:,i),'color','black','linewidth',2)
    hold on
end
for i = 1:size(bondary.LHW2_x,2)
    plot3(bondary.LHW2_x(:,i),bondary.LHW2_y(:,i),bondary.LHW2_z(:,i),'color','black','linewidth',2)
    hold on
end
nphi_limiter = atan(76/2350)/(pi/nphi);
for i = 1:nphi_limiter
    plot3(bondary.limiter_x(:,i),bondary.limiter_y(:,i),bondary.limiter_z(:,i),'color','black','linewidth',2)
end
% line(divertor_r,zeros(79,1),divertor_z)

xlabel('X', 'fontsize',18,'fontweight','bold');
ylabel('Y', 'fontsize',18,'fontweight','bold');
zlabel('Z', 'fontsize',18,'fontweight','bold');
axis equal

% title([efit.cas],'Fontsize',20);
if t(1,1) < 2*pi*(13/32)
    
    title(['left channel #',num2str(shotnum),' at ',num2str(tpoint),'s'],'Fontsize',20);
else 
    title(['right channel #',num2str(shotnum),' at ',num2str(tpoint),'s'],'Fontsize',20);
end
%---------------------plot EAST section ------------------------------

plot3(bondary.r(1:58,1),zeros(length(bondary.r(1:58,1)),1),bondary.z(1:58,1),'color','black','linewidth',4);hold on
plot3(bondary.r(58:71,1),zeros(length(bondary.r(58:71,1)),1),bondary.z(58:71,1),'color','c','linewidth',4)
plot3(bondary.r(71:90,1),zeros(length(bondary.r(71:90,1)),1),bondary.z(71:90,1),'color','black','linewidth',4)
% plot3(bondary.r(:,1),zeros(length(bondary.r(:,1)),1),bondary.z(:,1),'color','k','linewidth',4)
% bondary.r(56:71,i)=flip(bondary.LowerDivertor3(:,1));



figure(33)
% plot(Lc(1,:)',[1:length(Lc)],'linewidth',2);hold on
% plot(Lc(2,:)',[1:length(Lc)],'linewidth',2);hold on
plot(Lc(3,:)',[1:length(Lc)],'linewidth',2);hold on
if t(1,1) < 2*pi*(13/32)
%     plot(Lc(1,:)',[1:length(Lc)],'linewidth',2);hold on
    title(['left channel #',num2str(shotnum),' at ',num2str(tpoint),'s'],'Fontsize',20);
else
%     plot(Lc(2,:)',[1:length(Lc)],'linewidth',2);hold on
    title(['right channel #',num2str(shotnum),' at ',num2str(tpoint),'s'],'Fontsize',20);
end
% legend('left','right','total')

toc

% save()




















