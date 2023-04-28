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