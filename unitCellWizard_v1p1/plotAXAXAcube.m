function plotAXAXAcube(crystal,spatialHarmonics,thres)
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
    
    %plot unit cell
    %calculate the fill fraction
    ff=sum(UCB(:))/(crystal.unitcell.N1*crystal.unitcell.N2*crystal.unitcell.N3);

    %plot
    faceColSurf=[0.8 0.8 0.8];
    faceColCaps=[0.5 0.5 0.5];
    
    figure;
    UCS = smooth3(UCB,'gaussian');
    v=isosurface(YU,XU,ZU,UCS,0.5);
    h=patch(v,'FaceColor',faceColSurf);
    f=isocaps(YU,XU,ZU,UCS,0.5);
    h=patch(f,'FaceColor',faceColCaps);
%     view(10,10);
    axis('equal','tight');
    xlim([-crystal.unitcell.a/2 crystal.unitcell.a/2]);
    ylim([-crystal.unitcell.a/2 crystal.unitcell.a/2]);
    zlim([-crystal.unitcell.a/2 crystal.unitcell.a/2]);
%     camlight; lighting phong;
    title(['UNITCELL ff=' num2str(ff)]);
    grid('on'); box('on');
    
    UCBn=double(UCB);
    %plot faces
    figure
    subplot(3,3,1);hold all;
    pcolor(squeeze(XU(:,:,1)),squeeze(YU(:,:,1)),squeeze(GA(:,:,1)))
    shading flat;
    axis equal tight
    subplot(3,3,4);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),squeeze(GA(:,:,end)))
    shading flat;
    axis equal tight
    subplot(3,3,7);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),abs(squeeze(GA(:,:,end))-squeeze(GA(:,:,1))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,2);hold all;
    pcolor(squeeze(XU(:,1,:)),squeeze(ZU(:,1,:)),squeeze(GA(:,1,:)))
    shading flat;
    axis equal tight
    subplot(3,3,5);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),squeeze(GA(:,end,:)))
    shading flat;
    axis equal tight
    subplot(3,3,8);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),abs(squeeze(GA(:,end,:))-squeeze(GA(:,1,:))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,3);hold all;
    pcolor(squeeze(ZU(1,:,:)),squeeze(YU(1,:,:)),squeeze(GA(1,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,6);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),squeeze(GA(end,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,9);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),abs(squeeze(GA(end,:,:))-squeeze(GA(1,:,:))))
    shading flat; colorbar
    axis equal tight
    %plot faces
    figure
    subplot(3,3,1);hold all;
    pcolor(squeeze(XU(:,:,1)),squeeze(YU(:,:,1)),squeeze(UCBn(:,:,1)))
    shading flat;
    axis equal tight
    subplot(3,3,4);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),squeeze(UCBn(:,:,end)))
    shading flat;
    axis equal tight
    subplot(3,3,7);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),abs(squeeze(UCBn(:,:,end))-squeeze(UCBn(:,:,1))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,2);hold all;
    pcolor(squeeze(XU(:,1,:)),squeeze(ZU(:,1,:)),squeeze(UCBn(:,1,:)))
    shading flat;
    axis equal tight
    subplot(3,3,5);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),squeeze(UCBn(:,end,:)))
    shading flat;
    axis equal tight
    subplot(3,3,8);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),abs(squeeze(UCBn(:,end,:))-squeeze(UCBn(:,1,:))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,3);hold all;
    pcolor(squeeze(ZU(1,:,:)),squeeze(YU(1,:,:)),squeeze(UCBn(1,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,6);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),squeeze(UCBn(end,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,9);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),abs(squeeze(UCBn(end,:,:))-squeeze(UCBn(1,:,:))))
    shading flat; colorbar
    axis equal tight
    %plot faces
    figure
    subplot(3,3,1);hold all;
    pcolor(squeeze(XU(:,:,1)),squeeze(YU(:,:,1)),squeeze(UCS(:,:,1)))
    shading flat;
    axis equal tight
    subplot(3,3,4);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),squeeze(UCS(:,:,end)))
    shading flat;
    axis equal tight
    subplot(3,3,7);hold all;
    pcolor(squeeze(XU(:,:,end)),squeeze(YU(:,:,end)),abs(squeeze(UCBn(:,:,end))-squeeze(UCBn(:,:,1))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,2);hold all;
    pcolor(squeeze(XU(:,1,:)),squeeze(ZU(:,1,:)),squeeze(UCS(:,1,:)))
    shading flat;
    axis equal tight
    subplot(3,3,5);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),squeeze(UCS(:,end,:)))
    shading flat;
    axis equal tight
    subplot(3,3,8);hold all;
    pcolor(squeeze(XU(:,end,:)),squeeze(ZU(:,end,:)),abs(squeeze(UCBn(:,end,:))-squeeze(UCBn(:,1,:))))
    shading flat; colorbar
    axis equal tight
    subplot(3,3,3);hold all;
    pcolor(squeeze(ZU(1,:,:)),squeeze(YU(1,:,:)),squeeze(UCS(1,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,6);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),squeeze(UCS(end,:,:)))
    shading flat;
    axis equal tight
    subplot(3,3,9);hold all;
    pcolor(squeeze(ZU(end,:,:)),squeeze(YU(end,:,:)),abs(squeeze(UCS(end,:,:))-squeeze(UCS(1,:,:))))
    shading flat; colorbar
    axis equal tight
end