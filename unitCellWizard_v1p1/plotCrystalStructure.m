function plotCrystalStructure(ax,dat)
    cla(ax);
    hold(ax,'all');
    %plot primitve lattice vectors
    quiver3(ax,0,0,0,dat.t1(1),dat.t1(2),dat.t1(3),'Color','b');
    quiver3(ax,0,0,0,dat.t2(1),dat.t2(2),dat.t2(3),'Color','b');
    quiver3(ax,0,0,0,dat.t3(1),dat.t3(2),dat.t3(3),'Color','b');
    %plot reciprocal lattice vectors
    quiver3(ax,0,0,0,dat.T1(1),dat.T1(2),dat.T1(3),'Color','r');
    quiver3(ax,0,0,0,dat.T2(1),dat.T2(2),dat.T2(3),'Color','r');
    quiver3(ax,0,0,0,dat.T3(1),dat.T3(2),dat.T3(3),'Color','r');
    %plot path around IBZ
    plot3(ax,dat.KP(1,:),dat.KP(2,:),dat.KP(3,:),'k')
    KL=strsplit(dat.KL,' ');
    for ii=1:size(dat.KP,2)
        text(ax,dat.KP(1,ii),dat.KP(2,ii),dat.KP(3,ii),KL{ii});
    end
    view(ax,10,10)
    axis(ax,'equal','tight');
%             xlim(ax,[-app.a/2 app.a/2]);
%             ylim(ax,[-app.a/2 app.a/2]);
%             zlim(ax,[-app.a/2 app.a/2]);
    title(ax,'Crystal')
    grid(ax,'on'); box(ax,'on');
end