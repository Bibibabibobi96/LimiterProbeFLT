clear idx
for i = 1:size(flt.backward_r,2)
    idx_1 = find(isnan(flt.forward_x(:,i)), 1);
    if isempty(idx_1)
        idx(1,i) = flt.n * nphi;
    else
        idx(1,i) = idx_1;
    end
end

for i = 1:size(flt.backward_r,2)
    idx_1 = find(isnan(flt.backward_x(:,i)), 1);
    if isempty(idx_1)
        idx(2,i) = flt.n * nphi;
    else
        idx(2,i) = idx_1;
    end
end

idx(3,:) = idx(1,:) + idx(2,:);

%% Lc pcolor
load('G:\LimiterProbeFLT\EAST_bondary\EAST_new_divertor_limiter_3D.mat')
flt.rStart = linspace(min(flt.rStart_real), max(flt.rStart_real), flt.nr_mesh);
flt.zStart = linspace(min(flt.zStart_real), max(flt.zStart_real), flt.nz_mesh);

figure
t = tiledlayout(1,3);

ax(1) = nexttile(1);
% pcolor(flt.rStart,flt.zStart,reshape(idx(1,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
pcolor(flt.rStart,flt.zStart,reshape(flt.Lc(1,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
plot(bondary.r(1:23,1),bondary.z(1:23,1),'Color','k','LineWidth',2)
title('Forward','FontSize',14)
shading interp
colorbar
axis equal tight
ax(2) = nexttile(2);
% pcolor(flt.rStart,flt.zStart,reshape(idx(2,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
pcolor(flt.rStart,flt.zStart,reshape(flt.Lc(2,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
plot(bondary.r(1:23,1),bondary.z(1:23,1),'Color','k','LineWidth',2)
title('Backward','FontSize',14)
shading interp
colorbar
axis equal tight
ax(3) = nexttile(3);
% pcolor(flt.rStart,flt.zStart,reshape(idx(3,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
pcolor(flt.rStart,flt.zStart,reshape(flt.Lc(3,:) / nphi,flt.nz_mesh,flt.nr_mesh));hold on
plot(bondary.r(1:23,1),bondary.z(1:23,1),'Color','k','LineWidth',2)
title('Total','FontSize',14)
shading interp
colorbar
axis equal tight

t.TileSpacing = 'compact';
t.Padding = 'none';

sgtitle(['#',num2str(shot),' at ',num2str(tpoint),'s'],'FontSize',18)

%%

figure
t = tiledlayout(1,3);

nexttile
scatter(flt.backward_r(1:nphi:nphi*flt.n,:),flt.backward_z(1:nphi:nphi*flt.n,:),0.5);hold on
plot(bondary.r(:,1),bondary.z(:,1),'Color','k','LineWidth',2)  
plot(flt.backward_r(1,:),flt.backward_z(1,:),'o','Color','k','LineWidth',1)  
axis equal;axis tight;
title('Backward')

nexttile
scatter(flt.forward_r(1:nphi:nphi*flt.n,:),flt.forward_z(1:nphi:nphi*flt.n,:),0.5);hold on
plot(bondary.r(:,1),bondary.z(:,1),'Color','k','LineWidth',2) 
plot(flt.backward_r(1,:),flt.backward_z(1,:),'o','Color','k','LineWidth',1)  
axis equal;axis tight;
title('Forward')

nexttile
scatter(flt.forward_r(1:nphi:nphi*flt.n,:),flt.forward_z(1:nphi:nphi*flt.n,:),0.5);hold on
scatter(flt.backward_r(1:nphi:nphi*flt.n,:),flt.backward_z(1:nphi:nphi*flt.n,:),0.5);hold on
plot(bondary.r(:,1),bondary.z(:,1),'Color','k','LineWidth',2) 
plot(flt.backward_r(1,:),flt.backward_z(1,:),'o','Color','k','LineWidth',1)  
axis equal;axis tight;
title('Total')

t.TileSpacing = 'none';
t.Padding = 'none';
axis equal
axis tight




