% function [V1,F1,V2,F2]=generateLattice_v1(crystal,N1,N2,N3)
%     figure();
%     for n3=1:N3
%         for n2=1:N2
%             for n1=1:N1
%                 v=crystal.unitcell.surface.vertices;
%                 v(:,1)=v(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
%                 v(:,2)=v(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
%                 v(:,3)=v(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
% %                 h=patch('Faces',F,'Vertices',V);
%                 h=patch('Faces',crystal.unitcell.surface.faces,'Vertices',v,'FaceColor','g');
%                 set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
%                 v=crystal.unitcell.caps.vertices;
%                 v(:,1)=v(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
%                 v(:,2)=v(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
%                 v(:,3)=v(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
% %                 h=patch('Faces',F,'Vertices',V);
%                 h=patch('Faces',crystal.unitcell.caps.faces,'Vertices',v,'FaceColor','g');
%                 set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
%             end
%         end
%     end
% end
function [V1,F1,V2,F2]=generateLattice_v1(crystal,N1,N2,N3)
    figure();
    vss=crystal.unitcell.surface.vertices;
    vcs=crystal.unitcell.caps.vertices;
    vs=[];
    vc=[];
    for n3=-N3/2:N3/2
        for n2=-N2/2:N2/2
            for n1=-N1/2:N1/2
                
                vst(:,1)=vss(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
                vst(:,2)=vss(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
                vst(:,3)=vss(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
                vs=[vs;vst];
                vct(:,1)=vcs(:,1)+n1*crystal.t1(1)+n2*crystal.t2(1)+n3*crystal.t3(1);
                vct(:,2)=vcs(:,2)+n1*crystal.t1(2)+n2*crystal.t2(2)+n3*crystal.t3(2);
                vct(:,3)=vcs(:,3)+n1*crystal.t1(3)+n2*crystal.t2(3)+n3*crystal.t3(3);
                vc=[vc;vct];
            end
        end
    end
    h=patch('Faces',crystal.unitcell.surface.faces,'Vertices',vs,'FaceColor','g');
    set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
    h=patch('Faces',crystal.unitcell.caps.faces,'Vertices',vc,'FaceColor','g');
    set(h,'FaceColor',[0.5 0.7 0.5],'LineWidth',0.5);
    xlabel('x'); ylabel('y'); zlabel('z');
end