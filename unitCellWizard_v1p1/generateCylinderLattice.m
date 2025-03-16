function generateCylinderLattice(crystal,sH,radius,height,th)
    dx=crystal.unitcell.a/size(crystal.unitcell.XU,1);
    dy=crystal.unitcell.a/size(crystal.unitcell.XU,2);
    dz=crystal.unitcell.a/size(crystal.unitcell.XU,3);
    r=0:dx:radius; 
    phi=0:dx/radius:2*pi; 
    z=0:dz:height;
    [R,PHI,Z]=meshgrid(r,phi,z);
    X=R.*cos(PHI);
    Y=R.*sin(PHI);
    
    GA=zeros(size(X));
    
    sHnames=fieldnames(sH);
    for ii=1:length(sHnames)
        if isstruct(sH.(sHnames{ii}))
            h=sH.(sHnames{ii});
            Kx=h.t1;
            Ky=h.t2;
            Kz=h.t3;
            GAtemp=exp(1j*(Kx*X+Ky*Y+Kz*Z-h.phi*ones(size(Z))));
            GA=GA+h.A*GAtemp;
        end
    end
%     
    %clean up numerical noise
    GA=real(GA);
    GA=GA-min(GA(:));
    GA=GA/max(GA(:));
    
    ucs=GA>th;
%     ucs=smooth3(uc);
    figure;
    v=isosurface(X,Y,Z,ucs,0.5);
    h=patch(v,'FaceColor','g');
    f=isocaps(X,Y,Z,ucs,0.5);
    h=patch(f,'FaceColor','g');
    axis equal tight
    xlim(radius*[-1 1]); ylim(radius*[-1 1]); zlim([0 height]);
    grid on;
    box on;
    
    %ask if export to stl file
    answer=questdlg('Is lattice to be saved to stl?','Export lattice','Yes','No','No');
    if strcmp(answer,'Yes')
        [file,path] = uiputfile('*.stl','Export to STL');
    
        if and(~isnumeric(file),~isnumeric(path))
            [F,V]=isosurface(X,Y,Z,ucs,0.5);
            [FC,VC]=isocaps(X,Y,Z,ucs,0.5);
            %combine faces and vertices
            F=[F;FC+length(V(:,1))];
            V=[V;VC];
            %save to stl
            svlcad([path file],F,V);
        end
    end
end