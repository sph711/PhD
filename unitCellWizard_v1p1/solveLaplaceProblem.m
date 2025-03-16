function solveLaplaceProblem(uc,ga,th,nH,nL)
    UCB=ga>th;
    Nx=size(UCB,1); ax=uc.a/Nx;
    Ny=size(UCB,2); ay=uc.a/Ny;
    Nz=size(UCB,3); az=uc.a/Nz;

    %e field along z
    Cz=0;
    for ii=1:Nx
        for jj=1:Ny
            Ccol=0;
            for kk=1:Nz
                if UCB(ii,jj,kk)
                    Ctemp=ax*ay*nH^2/az;
                else
                    Ctemp=ax*ay*nL^2/az;
                end
                Ccol=Ccol+1/Ctemp;
            end
            Ccol=1/Ccol;
            Cz=Cz+Ccol;
        end
    end
    epz=Cz/uc.a;
    
    %e field along y
    Cy=0;
end