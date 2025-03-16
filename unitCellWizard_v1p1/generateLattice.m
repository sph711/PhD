function [V1,F1,V2,F2]=generateLattice(crystal,N1,N2,N3)
    figure();
    for n3=-N3/2:N3/2
        for n2=-N2/2:N2/2
            for n1=-N1/2:N1/2
                v=crystal.unitcell.surface.vertices;
                v(:,1)=v(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
                v(:,2)=v(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
                v(:,3)=v(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
%                 h=patch('Faces',F,'Vertices',V);
                h=patch('Faces',crystal.unitcell.surface.faces,'Vertices',v,'FaceColor','g');
                set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
                v=crystal.unitcell.caps.vertices;
                v(:,1)=v(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
                v(:,2)=v(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
                v(:,3)=v(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
%                 h=patch('Faces',F,'Vertices',V);
                h=patch('Faces',crystal.unitcell.caps.faces,'Vertices',v,'FaceColor','g');
                set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
            end
        end
    end
    axis equal tight
    view(50,15);
    camlight;
    xlabel('x'); ylabel('y'); zlabel('z');
    drawnow;
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