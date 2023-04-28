clear
clc

tic
% B_total;
shot=103950;
tpoint=5.5;


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
% 
% BR3D_total = Bfield.BR3D_total;
% BZ3D_total = Bfield.BZ3D_total;
% Bphi3D_total = Bfield.Bphi3D_total;

BR3D_total(:,:,nphi+1) = BR3D_total(:,:,2);
BZ3D_total(:,:,nphi+1) = BZ3D_total(:,:,2);
Bphi3D_total(:,:,nphi+1) = Bphi3D_total(:,:,2);

%% Field Line Tracing
clear flt
tic

flt.R_min = min(R);
flt.R_max = max(R);
flt.Z_min = min(Z);
flt.Z_max = max(Z);

% Define 3D bondary condition

% load('EAST_new_divertor.mat')
% load('EAST_new_divertor_limiter+0.01.mat')
% load('EAST_new_divertor_limiter+0.05_3D.mat')
% load('EAST_new_divertor_limiter+LHW+0.05_3D.mat')
load('G:\LimiterProbeFLT\EAST_bondary\EAST_new_divertor_limiter+0.01_3D.mat')
% load('G:\LimiterProbeFLT\EAST_bondary\EAST_square.mat')

flt.bondary_r = bondary.r;
flt.bondary_z = bondary.z;

% Define start points

% load ('G:\LimiterProbeFLT\limiterprobeposition.mat')
% load('G:\LimiterProbeFLT\gogogo137.mat')
% load ('FLTtest.mat')
% load ('limiterprobeposition_inner0.1.mat')
flt.n = 200;
flt.nr_mesh = 30;
flt.nz_mesh = 5;

flt.rStart = linspace(2.25, 2.35, flt.nr_mesh);
flt.zStart = linspace(-0.5, 0.5, flt.nz_mesh);
% flt.rStart = linspace(1.9, 2.35, flt.nr_mesh);
% flt.zStart = linspace(0, 0, flt.nz_mesh);

% flt.tStart=-t+2*pi*(13/32);
%correction to real EAST toroidal angel

% t=t+2*pi*(13/32);
[flt.rStart,flt.zStart] = meshgrid(flt.rStart,flt.zStart);
flt.rStart_real = reshape(flt.rStart,1, size(flt.rStart,1) * size(flt.rStart,2));
flt.zStart_real = reshape(flt.zStart,1, size(flt.zStart,1) * size(flt.zStart,2));
flt.tStart_real = 0.0335 * ones(1,flt.nr_mesh * flt.nz_mesh);
flt.nphiStart = ceil(flt.tStart_real(1)/(2*pi/(nphi-1)));


clear h
global h
h = 2 * pi / (nphi - 1);

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

% flt.Lc(2,:) = flt.Lc(1,:);

disp('Start Forward Tracing')
h = -h;

for i = 1 : flt.nr_mesh * flt.nz_mesh
    disp(['Start tracing from point',num2str(i)]);
    flt.rStart = flt.rStart_real(1,i);
    flt.zStart = flt.zStart_real(1,i);
    flt.tStart = flt.tStart_real(1,i);
    [flt.forward_r(:,i),flt.forward_z(:,i),flt.forward_t(:,i),flt.Lc(2,i)] = LPflt(flt,BR3D_total,BZ3D_total,Bphi3D_total);
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

Lc_pcolor

field_line_plot 




















