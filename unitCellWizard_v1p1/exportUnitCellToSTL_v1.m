function exportUnitCellToSTL_v1(uc,path,name,ucsComp)
    ucs = smooth3(uc.UCB);
    %make sure unitcell boundaries are continuous by assigning them the
    %mean of the values on both sides
    if ucsComp
        tmp=(ucs(end,:,:)+ucs(1,:,:))./2;
        ucs(end,:,:)=tmp; ucs(1,:,:)=tmp;
        tmp=(ucs(:,end,:)+ucs(:,1,:))./2;
        ucs(:,end,:)=tmp; ucs(:,1,:)=tmp;
        tmp=(ucs(:,:,end)+ucs(:,:,1))./2;
        ucs(:,:,end)=tmp; ucs(:,:,1)=tmp;
    end
    
    [F,V]=isosurface(uc.YU,uc.XU,uc.ZU,ucs,0.5);
    [FC,VC]=isocaps(uc.YU,uc.XU,uc.ZU,ucs,0.5);
    
    %combine faces and vertices
    F=[F;FC+length(V(:,1))];
    V=[V;VC];
    
    svlcad([path '\' name '.stl'],F,V);
end