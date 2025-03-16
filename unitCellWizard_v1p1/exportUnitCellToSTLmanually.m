function exportUnitCellToSTLmanually(uc)
    ucs= smooth3(uc.UCB);
    [F,V]=isosurface(uc.YU,uc.XU,uc.ZU,ucs,0.5);
    [FC,VC]=isocaps(uc.YU,uc.XU,uc.ZU,ucs,0.5);
    %combine faces and vertices
    F=[F;FC+length(V(:,1))];
    V=[V;VC];
    
    [file,path]=uiputfile('unitcell.stl');
    svlcad([path file],F,V);
end