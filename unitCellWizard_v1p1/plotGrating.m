function plotGrating(uc,ga,th)
    
figure
view(3);
hold all;
patch(surf2patch(squeeze(uc.XU(:,:,1)),squeeze(uc.YU(:,:,1)),squeeze(uc.ZU(:,:,1)),real(squeeze(ga(:,:,1)))))
patch(surf2patch(squeeze(uc.XU(:,:,end)),squeeze(uc.YU(:,:,end)),squeeze(uc.ZU(:,:,end)),real(squeeze(ga(:,:,end)))))
patch(surf2patch(squeeze(uc.XU(:,1,:)),squeeze(uc.YU(:,1,:)),squeeze(uc.ZU(:,1,:)),real(squeeze(ga(:,1,:)))))
patch(surf2patch(squeeze(uc.XU(:,end,:)),squeeze(uc.YU(:,end,:)),squeeze(uc.ZU(:,end,:)),real(squeeze(ga(:,end,:)))))
patch(surf2patch(squeeze(uc.XU(1,:,:)),squeeze(uc.YU(1,:,:)),squeeze(uc.ZU(1,:,:)),real(squeeze(ga(1,:,:)))))
patch(surf2patch(squeeze(uc.XU(end,:,:)),squeeze(uc.YU(end,:,:)),squeeze(uc.ZU(end,:,:)),real(squeeze(ga(end,:,:)))))
shading flat;colorbar; grid
caxis([0 1])
xlabel('x','Interpreter','latex');
ylabel('y','Interpreter','latex');
zlabel('z','Interpreter','latex');


%     %plot
%     figure;
%     UCS = smooth3(UCB);
%     v=isosurface(uc.YU,uc.XU,uc.ZU,UCS,0.5);
%     h=patch(v,'FaceColor',faceColSurf);
%     f=isocaps(uc.YU,uc.XU,uc.ZU,UCS,0.5);
%     h=patch(f,'FaceColor',faceColCaps);
%     view(30,25);
%     axis('equal','tight');
%     xlim([-uc.a/2 uc.a/2]);
%     ylim([-uc.a/2 uc.a/2]);
%     zlim([-uc.a/2 uc.a/2]);
%     set(gca,'TickLabelInterpreter','latex')
%     set(gca,'Fontsize',24);
%     set(gca,'XTickLabels',{'-a/2','0','a/2'})
%     set(gca,'YTickLabels',{'-a/2','0','a/2'})
%     set(gca,'ZTickLabels',{'-a/2','0','a/2'})
%     xlabel('x','Interpreter','latex');
%     ylabel('y','Interpreter','latex');
%     zlabel('z','Interpreter','latex');
% %     camlight; lighting('phong');
% %     title(['UNITCELL ff=' num2str(ff)]);
%     grid('on'); box('on');
%     
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end