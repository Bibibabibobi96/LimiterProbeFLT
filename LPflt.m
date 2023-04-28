function [ flt_r, flt_z, flt_t, Lc ] = LPflt(flt, BR3D_total, BZ3D_total, Bphi3D_total)
%flt struct should include:
% n: number of truns
% rStart zStart tStart :coordinates of start points
% nphi and h are global
% this program should be able to handel FLT with a series starting points
% BUT it only works with ONE SINGLE starting point in current version

% tic

clear x y p i
clear dBr dBz dBp dr dz
clear poincare_r poincare_z poincare_t
clear flt_r flt_z flt_t flt_x flt_y

bondary.r = flt.bondary_r;
bondary.z = flt.bondary_z;

% flt.n = 10; %number of truns
nr = size(BZ3D_total,1);
nz = size(BR3D_total,2);

% nphi=size(Bphi3D_total,3);
% h=2*pi/nphi;
global h nphi
R = linspace(flt.R_min,flt.R_max,nr);
Z = linspace(flt.Z_min,flt.Z_max,nz);

phi = linspace(0,2*pi,nphi);
% phi(nphi + 1) = phi(2)+2*pi;


% r = linspace(R_min,R_max,nr);
% z = linspace(Z_min,Z_max,nz);



r = flt.rStart;
z = flt.zStart;
t = flt.tStart;

R_max = flt.R_max;
R_min = flt.R_min;
Z_max = flt.Z_max;
Z_min = flt.Z_min;



%     disp(['Start tracing from point',num2str(n_point)]);

for ii = 1 : flt.n * nphi
    % for ii=1:289

    %         n_ind1 = mod(flt.nphiStart + sign(h) * ii -1 , nphi) + 1;
    %         n_ind2 = mod(flt.nphiStart + sign(h) * ii -1 , nphi) + 2;
    %  转一圈会少一个层
    %         if n_ind1 == 720
    %             n_ind1 = 1;
    %             n_ind2 = 3;
    %         end

    [~,index]=sort(abs(phi-mod(t(ii) + h/2,2*pi)));
    index = sort([index(1),index(2)]);

    n_ind1 = index(1);
    n_ind2 = index(2);



%         dBr(1) = interp2(R,Z,BR3D_total(:,:,n_ind1),r(ii),z(ii),'linear',nan);
%         dBz(1) = interp2(R,Z,BZ3D_total(:,:,n_ind1 ),r(ii),z(ii),'linear',nan);
%         dBp(1) = interp2(R,Z,Bphi3D_total(:,:,n_ind1),r(ii),z(ii),'linear',nan);
    dBr(1) = interp3_nearest_Liao(R,Z,phi,BR3D_total,r(ii),z(ii),t(ii),R_min,R_max,Z_min,Z_max);
    dBz(1) = interp3_nearest_Liao(R,Z,phi,BZ3D_total,r(ii),z(ii),t(ii),R_min,R_max,Z_min,Z_max);
    dBp(1) = interp3_nearest_Liao(R,Z,phi,Bphi3D_total,r(ii),z(ii),t(ii),R_min,R_max,Z_min,Z_max);
    dr(1) = dBr(1)./dBp(1);
    dz(1) = dBz(1)./dBp(1);


%         dBr(2) = interp3(R,Z,phi(n_ind1: n_ind2),BR3D_total(:,:,n_ind1: n_ind2),...
%             r(ii)+(h/2)*dr(1), z(ii)+(h/2)*dz(1), mod(t(ii)+h/2,2*pi),'linear',nan);
%         dBz(2) = interp3(R,Z,phi(n_ind1: n_ind2),BZ3D_total(:,:,n_ind1: n_ind2),...
%             r(ii)+(h/2)*dr(1), z(ii)+(h/2)*dz(1), mod(t(ii)+h/2,2*pi),'linear',nan);
%         dBp(2) = interp3(R,Z,phi(n_ind1: n_ind2),Bphi3D_total(:,:,n_ind1:n_ind2),...
%             r(ii)+(h/2)*dr(1), z(ii)+(h/2)*dz(1), mod(t(ii)+h/2,2*pi),'linear',nan);
    dBr(2) = interp3_nearest_Liao(R,Z,phi,BR3D_total,r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
    dBz(2) = interp3_nearest_Liao(R,Z,phi,BZ3D_total,r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
    dBp(2) = interp3_nearest_Liao(R,Z,phi,Bphi3D_total,r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
%     dBr(2) = interp2(R,Z,BR3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),'linear',nan);
%     dBz(2) = interp2(R,Z,BZ3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),'linear',nan);
%     dBp(2) = interp2(R,Z,Bphi3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(1),z(ii)+(h/2)*dz(1),'linear',nan);
    dr(2) = dBr(2)./dBp(2);
    dz(2) = dBz(2)./dBp(2);

%         dBr(3) = interp3(R,Z,phi(n_ind1:n_ind2),BR3D_total(:,:,n_ind1:n_ind2),...
%             r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),mod(t(ii)+h/2,2*pi),'linear',nan);
%         dBz(3) = interp3(R,Z,phi(n_ind1:n_ind2),BZ3D_total(:,:,n_ind1:n_ind2),...
%             r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),mod(t(ii)+h/2,2*pi),'linear',nan);
%         dBp(3) = interp3(R,Z,phi(n_ind1:n_ind2),Bphi3D_total(:,:,n_ind1:n_ind2),...
%             r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),mod(t(ii)+h/2,2*pi),'linear',nan);
    dBr(3) = interp3_nearest_Liao(R,Z,phi,BR3D_total,r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
    dBz(3) = interp3_nearest_Liao(R,Z,phi,BZ3D_total,r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
    dBp(3) = interp3_nearest_Liao(R,Z,phi,Bphi3D_total,r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),t(ii)+h/2,R_min,R_max,Z_min,Z_max);
%     dBr(3) = interp2(R,Z,BR3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),'linear',nan);
%     dBz(3) = interp2(R,Z,BZ3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),'linear',nan);
%     dBp(3) = interp2(R,Z,Bphi3D_total(:,:,n_ind1),r(ii)+(h/2)*dr(2),z(ii)+(h/2)*dz(2),'linear',nan);
    dr(3) = dBr(3)./dBp(3);
    dz(3) = dBz(3)./dBp(3);

    dBr(4) = interp3_nearest_Liao(R,Z,phi,BR3D_total,r(ii)+h*dr(3),z(ii)+h*dz(3),t(ii)+h,R_min,R_max,Z_min,Z_max);
    dBz(4) = interp3_nearest_Liao(R,Z,phi,BZ3D_total,r(ii)+h*dr(3),z(ii)+h*dz(3),t(ii)+h,R_min,R_max,Z_min,Z_max);
    dBp(4) = interp3_nearest_Liao(R,Z,phi,Bphi3D_total,r(ii)+h*dr(3),z(ii)+h*dz(3),t(ii)+h,R_min,R_max,Z_min,Z_max);
%     dBr(4) = interp2(R,Z,BR3D_total(:,:,n_ind2),r(ii)+h*dr(3),z(ii)+h*dz(3),'linear',nan);
%     dBz(4) = interp2(R,Z,BZ3D_total(:,:,n_ind2),r(ii)+h*dr(3),z(ii)+h*dz(3),'linear',nan);
%     dBp(4) = interp2(R,Z,Bphi3D_total(:,:,n_ind2),r(ii)+h*dr(3),z(ii)+h*dz(3),'linear',nan);
    dr(4) = dBr(4)./dBp(4);
    dz(4) = dBz(4)./dBp(4);


    r(ii+1) = r(ii) + (h/6)*(dr(1)+2*dr(2)+2*dr(3)+dr(4));
    z(ii+1) = z(ii) + (h/6)*(dz(1)+2*dz(2)+2*dz(3)+dz(4));
    t(ii+1) = t(ii) + h;

    % clear dBr dBz dBp
    %%%%---find out wether current step is still in the EAST or not----

    if sum(isnan(dr)) == 0
        if h>0 && mod(ii + flt.nphiStart,nphi) == 0
            in_or_out = inpolygon(r(ii+1),z(ii+1),bondary.r(:,nphi+1),bondary.z(:,nphi+1));
        elseif h<0 && mod(-ii + flt.nphiStart,nphi) == 0
            in_or_out = inpolygon(r(ii+1),z(ii+1),bondary.r(:,nphi+1),bondary.z(:,nphi+1));
        end

        if h>0 && mod(ii + flt.nphiStart,nphi) ~= 0
            in_or_out = inpolygon(r(ii+1),z(ii+1),bondary.r(:,mod(ii + flt.nphiStart,nphi)),bondary.z(:,mod(ii + flt.nphiStart,nphi)));
        elseif h<0 && mod(-ii + flt.nphiStart,nphi) ~= 0
            in_or_out = inpolygon(r(ii+1),z(ii+1),bondary.r(:,mod(-ii + flt.nphiStart,nphi)),bondary.z(:,mod(-ii + flt.nphiStart,nphi)));
        end

    else
        in_or_out = 0;
    end


    if in_or_out == 0

        Lc(1) = ii;
        disp(['End tracing at step ',num2str(ii)]);
        r(ii+1:flt.n*nphi+1) = NaN;
        z(ii+1:flt.n*nphi+1) = NaN;
        t(ii+1:flt.n*nphi+1) = NaN;
        break

    elseif ii == flt.n*nphi
        Lc(1)=ii;
    end

    %%%%%--------------end----------------
    if rem(ii,nphi)==0
        t(ii+1)=0;
        t(ii+1)=t(ii)+h;
    end

end



flt_r = r';
flt_t = t';
flt_z = z';

end
% clear r t z

% toc





