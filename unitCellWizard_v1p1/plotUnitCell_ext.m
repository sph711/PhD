function plotUnitCell_ext(uc,ga,th)
    faceColSurf=[0.8 0.8 0.8]
    faceColCaps=[0.5 0.5 0.5]
    %generate binary unit cell
    UCB=ga>th;

    %calculate the fill fraction
    ff=sum(UCB(:))/(uc.N1*uc.N2*uc.N3);

    %plot
    figure;
    UCS = smooth3(UCB);
    v=isosurface(uc.YU,uc.XU,uc.ZU,UCS,0.5);
    h=patch(v,'FaceColor',faceColSurf);
    f=isocaps(uc.YU,uc.XU,uc.ZU,UCS,0.5);
    h=patch(f,'FaceColor',faceColCaps);
    view(30,25);
    axis('equal','tight');
    xlim([-uc.a/2 uc.a/2]);
    ylim([-uc.a/2 uc.a/2]);
    zlim([-uc.a/2 uc.a/2]);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'Fontsize',24);
    set(gca,'XTickLabels',{'$-a/2$','$0$','$a/2$'})
    set(gca,'YTickLabels',{'$-a/2$','$0$','$a/2$'})
    set(gca,'ZTickLabels',{'$-a/2$','$0$','$a/2$'})
    xlabel('x','Interpreter','latex');
    ylabel('y','Interpreter','latex');
    zlabel('z','Interpreter','latex');
%     camlight; lighting('phong');
%     title(['UNITCELL ff=' num2str(ff)]);
    grid('on'); box('on');
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end