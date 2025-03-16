function [UCB,v,f]=plotUnitCell(ax,uc,ga,th)
%     ga=ifftshift(ga);
    %generate binary unit cell
    UCB=ga>th;

    %calculate the fill fraction
    ff=sum(UCB(:))/(uc.N1*uc.N2*uc.N3);

    %plot
    cla(ax);
%     UCS = smooth3(UCB);
    UCS=UCB;
    v=isosurface(uc.YU,uc.XU,uc.ZU,UCS,0.5);
    h=patch(ax,v,'FaceColor','g');
    f=isocaps(uc.YU,uc.XU,uc.ZU,UCS,0.5);
    h=patch(ax,f,'FaceColor','g');
    view(ax,10,10);
    axis(ax,'equal','tight');
    xlim(ax,[-uc.a/2 uc.a/2]);
    ylim(ax,[-uc.a/2 uc.a/2]);
    zlim(ax,[-uc.a/2 uc.a/2]);
    camlight(ax); lighting(ax,'phong');
    title(ax,['UNITCELL ff=' num2str(ff)]);
    grid(ax,'on'); box(ax,'on');
end