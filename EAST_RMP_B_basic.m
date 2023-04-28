% coil_geometry_circle2
% clear
clear all
% load('high_m_coordinates.mat')
load('G:\LimiterProbeFLT\RMPxyz.mat')
load(['G:\ERGOS_HCFs\dig\notebooks\103950_5000.0.mat'])
% clear
clear BZ3D
clear BR3D
clear Bphi3D
clear Br_3D Bz_3D Bphi_3D
% load('RMPxyz')
% load('param.dat')

% global IU IL
global nphi
% 

for i = 1:16
    coils_x{i} = RMPxyz{i}(1,:);
    coils_y{i} = RMPxyz{i}(2,:);
    coils_z{i} = RMPxyz{i}(3,:);
end
I=1; 
% IU=I
% IL=I
        
%   nR=param(1);
%   nZ=param(2);
%   nphi=param(3);
%   Rmin=param(4);
%   Rmax=param(5);
%   Zmin=param(6);
%   Zmax=param(7);
  

k=0;

% Br=0;Bz=0;Bphi=0;
%   nR=180;
%   nZ=257;
%   nphi=360;
  nR=129;
  nZ=129;
  nphi=720;


  for i = 1:16
      coils_x{i} = permute(coils_x{i},[2 1]);
      coils_y{i} = permute(coils_y{i},[2 1]);
      coils_z{i} = permute(coils_z{i},[2 1]);
  end

% R = linspace(Rmin,Rmax,nR);
% Z = linspace(Zmin,Zmax,nZ);
% phi = linspace(0,2*pi,nphi);
% [R_3D Z_3D phi_3D] = meshgrid(R,Z,phi);

[R_3D Z_3D ] = meshgrid(R,Z);

% I = gpuArray(I);
k=0;  
n=1;

for i = 1:16
    BR3D{i}=zeros(nZ,nR,nphi);
    BZ3D{i}=zeros(nZ,nR,nphi);
    Bphi3D{i}=zeros(nZ,nR,nphi);
end

h=waitbar(0,'please wait');

clear Br Bz Bphi
% Br{1} = 0; 
% Bz{1} = 0;
% Bphi{1} = 0;

for phi = 0:2*pi/(nphi-1):2*pi
    k=k+1;
    
    
    for j = 1:16
        for i = 1:n:length(coils_x{5})-n
            
            
            [Br0,Bz0,Bphi0]=finiteelementC(I,(coils_x{j}(i)),coils_y{j}(i),coils_z{j}(i),...
                coils_x{j}(i+n),coils_y{j}(i+n),coils_z{j}(i+n),R_3D,Z_3D,phi);
            if i == 1
                Br{j} = Br0;
                Bz{j} = Bz0;
                Bphi{j} = Bphi0;
                
            else
                Br{j} = Br{j}+Br0;
                Bz{j} = Bz{j}+Bz0;
                Bphi{j} = Bphi{j}+Bphi0;
               
            end
            
        end
        
        
        
    end
%     [i,j,k]
    for i = 1:16
        BR3D_RMP{i}(:,:,k) = Br{i};
        BZ3D_RMP{i}(:,:,k) = Bz{i};
        Bphi3D_RMP{i}(:,:,k) = Bphi{i};
    end
    
    
    str=['3D magnetic field calculating...',num2str(round(k/nphi*10000)/100),'%'];
    waitbar(k/nphi,h,str)
end
delete(h)



save (['G:\LimiterProbeFLT\B_field\', 'RMP_basic.mat'],'BR3D_RMP','BZ3D_RMP','Bphi3D_RMP')






