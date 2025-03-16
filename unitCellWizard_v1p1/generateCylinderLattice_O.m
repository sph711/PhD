function generateCylinderLattice_O(crystal,sH,radius,height,th)
    dx=crystal.unitcell.a/size(crystal.unitcell.XU,1);
    dy=crystal.unitcell.a/size(crystal.unitcell.XU,2);
    dz=crystal.unitcell.a/size(crystal.unitcell.XU,3);
    
    x=-radius:dx:radius;
    y=-radius:dy:radius;
    z=-height/2:dz:height/2;
    [X,Y,Z]=meshgrid(x,y,z);
%     
%     X(sqrt(X.^2+Y.^2)>radius)=NaN;
%     Y(sqrt(X.^2+Y.^2)>radius)=NaN;
%     Z(sqrt(X.^2+Y.^2)>radius)=NaN;

%     [X,Y,Z]=cylinder([1 1]*radius);

%transfrom all points outside the cylinder to the cylinder shell
    ind=find(sqrt(X.^2+Y.^2)>radius);
    for ii=1:length(ind)
        phi=atan2(Y(ind(ii)),X(ind(ii)));
        X(ind(ii))=radius*cos(phi);
        Y(ind(ii))=radius*sin(phi);
    end
    
    %rotate mesh
    theta=45*pi/180;
    phi=0*pi/180;
    gamma=45*pi/180;
    
    Xrot=zeros(size(X));
    Yrot=zeros(size(Y));
    Zrot=zeros(size(Z));
%    
    c1=cos(theta);
    c2=cos(phi);
    c3=cos(gamma);
    s1=sin(theta);
    s2=sin(phi);
    s3=sin(gamma);
    
    R_zyx=[c1*c2 c1*s2*s3-c3*s1 s1*s3+c1*c3*s2;
    c2*s1 c1*c3+s1*s2*s3 c3*s1*s2-c1*s3;
    -s2 c2*s3 c2*c3];

    temp=[X(:),Y(:),Z(:)]*R_zyx.' ;
      sz=size(X);
    Xrot=reshape(temp(:,1),sz);
    Yrot=reshape(temp(:,2),sz);
    Zrot=reshape(temp(:,3),sz);
    
    X=Xrot;
    Y=Yrot;
    Z=Zrot;
%     X(ind)=0;
%     Y(ind)=0;
%     Z(ind)=0;
%     
%     [R,PHI,Z]=meshgrid(r,phi,z);
%     X=R.*cos(PHI);
%     Y=R.*sin(PHI);
    
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
    ucs=smooth3(ucs);
    
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
            
            %clean up NaN values
            for ii=1:size(V,1)
                vToDel=[];
                if any(isnan(V(ii,:)))                    
                    %delete all faces with vertice
                    [i,~]=find(F==ii);
                    if ~isempty(i)
                        F(i,:)=[];
                    end
%                     while ~isempty(i)
%                         F(i(1),:)=[];
%                         [i,~]=find(F==ii);
%                     end
                end
            end
            %save to stl
            svlcad([path file],F,V);
        end
    else
            figure;
        v=isosurface(X,Y,Z,ucs,0.5);
        h=patch(v,'FaceColor','g');
        f=isocaps(X,Y,Z,ucs,0.5);
        h=patch(f,'FaceColor','g');
        axis equal tight
    %     xlim(radius*[-1 1]); ylim(radius*[-1 1]); zlim([0 height]);
        grid on;
        box on;
        xlabel('x');ylabel('y');zlabel('z');

    end
end