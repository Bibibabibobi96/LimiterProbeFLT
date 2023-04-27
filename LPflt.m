function [ flt_r, flt_z, flt_t, Lc ] = LPflt(flt, BR3D_total, BZ3D_total, Bphi3D_total)
%flt struct should include:
% n: number of truns
% R_min, R_amx
% Z_min, Z_max
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
R_max = flt.R_max;
R_min = flt.R_min;
Z_max = flt.Z_max;
Z_min = flt.Z_min;

phi = linspace(0,2*pi,nphi);

% r = linspace(R_min,R_max,nr);
% z = linspace(Z_min,Z_max,nz);



n_point=0;

r = flt.rStart;
z = flt.zStart;
t = flt.tStart;



for i = 1:size(r,2)

    jj = n_point+1;
    n_point = n_point+1;

%     disp(['Start tracing from point',num2str(n_point)]);

    for ii = 1 : flt.n * nphi
        % for ii=1:289

%         if isnan(r(ii,jj)) ~= 1 && isnan(z(ii,jj)) ~= 1
        if ~isnan(r(ii,jj)) && ~isnan(z(ii,jj))
            dBr(1,i) = interp3_linear_Liao(R,Z,phi,BR3D_total,r(ii,jj),z(ii,jj),t(ii,jj),R_min,R_max,Z_min,Z_max);
            dBz(1,i) = interp3_linear_Liao(R,Z,phi,BZ3D_total,r(ii,jj),z(ii,jj),t(ii,jj),R_min,R_max,Z_min,Z_max);
            dBp(1,i) = interp3_linear_Liao(R,Z,phi,Bphi3D_total,r(ii,jj),z(ii,jj),t(ii,jj),R_min,R_max,Z_min,Z_max);
            dr(1,i) = dBr(1,i)./dBp(1,i);
            dz(1,i) = dBz(1,i)./dBp(1,i);
        else
            dr(1,i) = NaN;  dz(1,i) = NaN;
        end

%         if isnan(dr(1)) ~= 1 && isnan(dz(1)) ~= 1
        if ~isnan(dr(1,i)) && ~isnan(dz(1,i))
            dBr(2,i) = interp3_linear_Liao(R,Z,phi,BR3D_total,r(ii,jj)+(h/2)*dr(1),z(ii,jj)+(h/2)*dz(1),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dBz(2,i) = interp3_linear_Liao(R,Z,phi,BZ3D_total,r(ii,jj)+(h/2)*dr(1),z(ii,jj)+(h/2)*dz(1),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dBp(2,i) = interp3_linear_Liao(R,Z,phi,Bphi3D_total,r(ii,jj)+(h/2)*dr(1),z(ii,jj)+(h/2)*dz(1),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dr(2,i) = dBr(2,i)./dBp(2,i);
            dz(2,i) = dBz(2,i)./dBp(2,i);
        else
            dr(2,i) = NaN;  dz(2,i) = NaN;
        end

%         if isnan(dr(2)) ~= 1 && isnan(dz(2)) ~= 1
        if ~isnan(dr(2,i)) && ~isnan(dz(2,i))
            dBr(3,i) = interp3_linear_Liao(R,Z,phi,BR3D_total,r(ii,jj)+(h/2)*dr(2),z(ii,jj)+(h/2)*dz(2),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dBz(3,i) = interp3_linear_Liao(R,Z,phi,BZ3D_total,r(ii,jj)+(h/2)*dr(2),z(ii,jj)+(h/2)*dz(2),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dBp(3,i) = interp3_linear_Liao(R,Z,phi,Bphi3D_total,r(ii,jj)+(h/2)*dr(2),z(ii,jj)+(h/2)*dz(2),t(ii,jj)+h/2,R_min,R_max,Z_min,Z_max);
            dr(3,i) = dBr(3,i)./dBp(3,i);
            dz(3,i) = dBz(3,i)./dBp(3,i);
        else
            dr(3,i) = NaN;  dz(3,i) = NaN;
        end

%         if isnan(dr(3)) ~= 1 && isnan(dz(3)) ~= 1
        if ~isnan(dr(3,i)) && ~isnan(dz(3,i))
            dBr(4,i) = interp3_linear_Liao(R,Z,phi,BR3D_total,r(ii,jj)+h*dr(3),z(ii,jj)+h*dz(3),t(ii,jj)+h,R_min,R_max,Z_min,Z_max);
            dBz(4,i) = interp3_linear_Liao(R,Z,phi,BZ3D_total,r(ii,jj)+h*dr(3),z(ii,jj)+h*dz(3),t(ii,jj)+h,R_min,R_max,Z_min,Z_max);
            dBp(4,i) = interp3_linear_Liao(R,Z,phi,Bphi3D_total,r(ii,jj)+h*dr(3),z(ii,jj)+h*dz(3),t(ii,jj)+h,R_min,R_max,Z_min,Z_max);
            dr(4,i) = dBr(4,i)./dBp(4,i);
            dz(4,i) = dBz(4,i)./dBp(4,i);
        else
            dr(4) = NaN;  dz(4) = NaN;
        end


%         r(ii+1,:)=r(ii,:)+(h/6)*(dr(1)+2*dr(2)+2*dr(3)+dr(4));
%         z(ii+1,:)=z(ii,:)+(h/6)*(dz(1)+2*dz(2)+2*dz(3)+dz(4));
        r(ii+1,jj)=r(ii,jj)+(h/6)*(dr(1)+2*dr(2)+2*dr(3)+dr(4));
        z(ii+1,jj)=z(ii,jj)+(h/6)*(dz(1)+2*dz(2)+2*dz(3)+dz(4));
        t(ii+1,jj)=t(ii,jj)+h;

        % clear dBr dBz dBp
        %%%%---find out wether current step is still in the EAST or not----

        if sum(isnan(dr))==0
            if h>0 && mod(ii + flt.nphiStart,nphi)==0
%                 [in_or_out] = in_or_out_area(bondary.r(:,nphi+1),bondary.z(:,nphi+1),r(ii+1,jj),z(ii+1,jj));
                [in_or_out] = inpolygon(r(ii+1,jj),z(ii+1,jj),bondary.r(:,nphi+1),bondary.z(:,nphi+1));
            elseif h<0 && mod(-ii + flt.nphiStart,nphi)==0
%                 [in_or_out] = in_or_out_area(bondary.r(:,nphi+1),bondary.z(:,nphi+1),r(ii+1,jj),z(ii+1,jj));
                [in_or_out] = inpolygon(r(ii+1,jj),z(ii+1,jj),bondary.r(:,nphi+1),bondary.z(:,nphi+1));
            end

            if h>0 && mod(ii + flt.nphiStart,nphi)~=0
%                 [in_or_out] = in_or_out_area(bondary.r(:,mod(ii + flt.nphiStart,nphi)),bondary.z(:,mod(ii + flt.nphiStart,nphi)),r(ii+1,jj),z(ii+1,jj));
                [in_or_out] = inpolygon(r(ii+1,jj),z(ii+1,jj),bondary.r(:,mod(ii + flt.nphiStart,nphi)),bondary.z(:,mod(ii + flt.nphiStart,nphi)));
            elseif h<0 && mod(-ii + flt.nphiStart,nphi)~=0
%                 [in_or_out] = in_or_out_area(bondary.r(:,mod(-ii + flt.nphiStart,nphi)),bondary.z(:,mod(-ii + flt.nphiStart,nphi)),r(ii+1,jj),z(ii+1,jj));
                [in_or_out] = inpolygon(r(ii+1,jj),z(ii+1,jj),bondary.r(:,mod(-ii + flt.nphiStart,nphi)),bondary.z(:,mod(-ii + flt.nphiStart,nphi)));
            end

        else
            in_or_out = 0;
        end



        if in_or_out == 0

            Lc(1,i)=ii;
            disp(['End tracing at step ',num2str(ii)]);
            r(ii+1:flt.n*nphi+1,jj) = NaN;
            z(ii+1:flt.n*nphi+1,jj) = NaN;
            t(ii+1:flt.n*nphi+1,jj) = NaN;
            break
            return

            %              if n_point == 1
            %               y = input('qingshuru')
            %              end
        elseif ii == flt.n*nphi
            Lc(1,i)=ii;
        end

        %%%%%--------------end----------------
        if rem(ii,nphi)==0
            t(ii+1,1:jj)=0;
            t(ii+1,1:jj)=t(ii,1:jj)+h;
        end

    end





    % poincare_r(1:ii+1,i)=r;
    % poincare_t(1:ii+1,i)=t;
    % poincare_z(1:ii+1,i)=z;

    flt_r(1:size(r,1),i)=r(:,jj);
    flt_t(1:size(t,1),i)=t(:,jj);
    flt_z(1:size(z,1),i)=z(:,jj);

end
% clear r t z

% toc





