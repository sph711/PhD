function uc=getAXAXAcube(crystal,spatialHarmonics,thres)
    GA=zeros(crystal.unitcell.N1,crystal.unitcell.N2,crystal.unitcell.N3);
    
    t1=crystal.unitcell.a*[1; 0; 0];
    t2=crystal.unitcell.a*[0; 1; 0];
    t3=crystal.unitcell.a*[0; 0; 1];
    
    p=linspace(-0.5,0.5,crystal.unitcell.N1);
    q=linspace(-0.5,0.5,crystal.unitcell.N2);
    r=linspace(-0.5,0.5,crystal.unitcell.N3);
    [Q,P,R]=meshgrid(q,p,r);
    XU=P*t1(1)+Q*t2(1)+R*t3(1);
    YU=P*t1(2)+Q*t2(2)+R*t3(2);
    ZU=P*t1(3)+Q*t2(3)+R*t3(3);
   
    sHnames=fieldnames(spatialHarmonics);
    for ii=1:length(sHnames)
        if isstruct(spatialHarmonics.(sHnames{ii}))
            h=spatialHarmonics.(sHnames{ii});
            Kx=h.t1;
            Ky=h.t2;
            Kz=h.t3;
            GAtemp=exp(1j*(Kx*XU+Ky*YU+Kz*ZU-h.phi*ones(size(ZU))));
            GA=GA+h.A*GAtemp;
        end
    end

    %clean up numerical noise
    GA=real(GA);
    GA=GA-min(GA(:));
    GA=GA/max(GA(:));
    
    %generate binary unit cell
    UCB=GA>thres;
    
    %make sure unitcell boundaries are equal
    UCB(end,:,:)=UCB(1,:,:);
    UCB(:,end,:)=UCB(:,1,:);
    UCB(:,:,end)=UCB(:,:,1);
    
    %result struct
    uc.UCB=UCB;
    uc.XU=XU;
    uc.YU=YU;
    uc.ZU=ZU;
end