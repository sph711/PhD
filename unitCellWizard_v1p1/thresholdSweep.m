function thresholdSweep(nth,cr,ga)
    th=linspace(0,1,nth);
    ff=zeros(1,nth);
    I1=ones(1,nth);
    I2=ones(1,nth);
    I3=ones(1,nth);
    for ii=1:nth
        UCB=ga>th(ii);
        ff(ii)=sum(UCB(:))/(cr.unitcell.N1*cr.unitcell.N2*cr.unitcell.N3);
        %detect if unit cell boundaries are empty
        s=isocaps(cr.unitcell.YU,cr.unitcell.XU,cr.unitcell.ZU,smooth3(UCB),0.5);
        if ~isempty(s.faces)
            N=zeros(size(s.faces));
            for jj=1:size(s.faces,1)
                f=s.faces(jj,:);
                v1=s.vertices(f(1),:);
                v2=s.vertices(f(2),:);
                v3=s.vertices(f(3),:);
                %compute unit surface normal
                n=cross((v2-v1),(v3-v1))/(norm(cross((v2-v1),(v3-v1))));
                N(jj,:)=n;
            end
            Nu=unique(N,'rows')
            if size(Nu,1)<6
                nt1=cr.t1/norm(cr.t1);
                if ~ismember(nt1',Nu,'rows')
                    I1(ii)=0;
                end
                nt2=cr.t2/norm(cr.t2);
                if ~ismember(nt2',Nu,'rows')
                    I2(ii)=0;
                end
                nt3=cr.t3/norm(cr.t3);
                if ~ismember(nt3',Nu,'rows')
                    I3(ii)=0;
                end
            end
        else
            I1(ii)=0;
            I2(ii)=0;
            I3(ii)=0;
        end
        %compute if the unit cell boundary is empty by calculating all face
        %normals see if all unique face normals correspond to primiteve
        %lattice vectors
    end
    figure; hold all
    plot(th,ff,'LineWidth',2);
    plot(th,I1,'LineWidth',2);
    plot(th,I2,'LineWidth',2);
    plot(th,I3,'LineWidth',2);
%     plot([th(I) th(I)],[0 1],'k--');
    grid on; box on;
    xlabel('Threshold'); ylabel('Fill Fraction');
    
end