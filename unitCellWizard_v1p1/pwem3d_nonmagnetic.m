function [BETA,WN]=pwem3d_nonmagnetic(crystal,pwemParam,sh,th)
    figure();
    NSOL=pwemParam.Nsol;
    wn_max=2;
    %% calculate beta axis resulution
    L=0;
    NKP=length(crystal.KP(1,:));
    for nkp=1:NKP-1
        L=L+norm(crystal.KP(:,nkp+1)-crystal.KP(:,nkp));
    end
    res=L/pwemParam.Nbeta;

    %% build bloch vector list
    BETA = crystal.KP(:,1);
    KT=1;
    for nkp=1:NKP-1
        %get endpoints of path
        kp1=crystal.KP(:,nkp);
        kp2=crystal.KP(:,nkp+1);
        %calculate number of points
        Lp=norm(kp2-kp1);
        NB = round(Lp/res);
        %calculate beta points from kp1 to kp2
        bx = kp1(1)+(kp2(1)-kp1(1))*[1:NB]/NB;
        by = kp1(2)+(kp2(2)-kp1(2))*[1:NB]/NB;
        bz = kp1(3)+(kp2(3)-kp1(3))*[1:NB]/NB;
        %append points to beta
        BETA=[BETA , [bx;by;bz]];
        KT(nkp+1)=length(BETA(1,:));
    end
    crystal.KL=strsplit(crystal.KL,' ');
    NBETA=length(BETA(1,:));

%     %show path in IBZ
%     figure
%     line(BETA(1,:),BETA(2,:),BETA(3,:),'LineWidth',2);
%     axis equal tight;
%     view(110,20)
    
    %% build the unit cell on fine qrid
    %construct an oblique meshgrid
    p=linspace(-0.5,0.5,pwemParam.N1);
    q=linspace(-0.5,0.5,pwemParam.N2);
    r=linspace(-0.5,0.5,pwemParam.N3);
    [Q,P,R]=meshgrid(q,p,r);
    XU=P*crystal.t1(1)+Q*crystal.t2(1)+R*crystal.t3(1);
    YU=P*crystal.t1(2)+Q*crystal.t2(2)+R*crystal.t3(2);
    ZU=P*crystal.t1(3)+Q*crystal.t2(3)+R*crystal.t3(3);
    
    %allocate space for 3d grating
    GA=zeros(pwemParam.N1,pwemParam.N2,pwemParam.N3);
    
    %compute 3d grating
    sHnames=fieldnames(sh);
    for ii=1:length(sHnames)
        if isstruct(sh.(sHnames{ii}))
            h=sh.(sHnames{ii});
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
    UC=GA>th;
    
    %compute fill fraction
    FF=sum(sum(sum(UC)))/prod(size(UC));
    
    %compute effective unit cell permittivity according to maxwell garnett
    epEffMG=pwemParam.nH^2*(2*FF*(pwemParam.nL^2-pwemParam.nH^2)+pwemParam.nL^2+2*pwemParam.nH^2)/(2*pwemParam.nH^2+pwemParam.nL^2-FF*(pwemParam.nL^2-pwemParam.nH^2));
    
    %convert unit cell to real materials
    ER=pwemParam.nL^2 +(pwemParam.nH^2-pwemParam.nL^2)*UC;
    ERMG=ones(size(ER))*epEffMG;
    %% compute photonic bands
    NP=pwemParam.NP; NQ=pwemParam.NQ; NR=pwemParam.NR;
    
    %build convolution matrix
    ERC=convmat(ER,NP,NQ,NR);
    ERMGC=convmat(ERMG,NP,NQ,NR);
    
    %initialize data records
    NH=NQ*NP*NR;
    WN=zeros(NSOL,NBETA);

    %calculate fourier expansion coefficiens
    p=[-floor(NP/2):floor(NP/2)];
    q=[-floor(NQ/2):floor(NQ/2)];
    r=[-floor(NR/2):floor(NR/2)];
    [Q,P,R]=meshgrid(q,p,r);
    
    %main loop - iterate over bloch wave vector
    for nb=1:NBETA

        %calculate wave vector expansion
        Kx=BETA(1,nb)-P*crystal.T1(1)-Q*crystal.T2(1)-R*crystal.T3(1);
        Ky=BETA(2,nb)-P*crystal.T1(2)-Q*crystal.T2(2)-R*crystal.T3(2);
        Kz=BETA(3,nb)-P*crystal.T1(3)-Q*crystal.T2(3)-R*crystal.T3(3);
        K=sqrt(abs(Kx).^2+abs(Ky).^2+abs(Kz).^2);

        %calculate perpendicular polarization vectors
        p1x=speye(NH,NH); p2x=speye(NH,NH);
        p1y=speye(NH,NH); p2y=speye(NH,NH);
        p1z=speye(NH,NH); p2z=speye(NH,NH);
        for m=1:NH
            k=[Kx(m);Ky(m);Kz(m)];
            if norm(k)<1e-10
                p1=[1;0;0];
                p2=[0;1;0];
            else
                v=[4*k(2); 2*k(3); 3*k(1)];
                p1=cross(k,v);
                p1=p1/norm(p1);
                p2=cross(k,p1);
                p2=p2/norm(p2);
            end
            p1x(m,m)=p1(1);
            p1y(m,m)=p1(2);
            p1z(m,m)=p1(3);
            p2x(m,m)=p2(1);
            p2y(m,m)=p2(2);
            p2z(m,m)=p2(3);
        end

        % build the eigenvalue problem
        K=diag(sparse(K(:)));
        
        A22=+K*p2x/ERC*K*p2x +K*p2y/ERC*K*p2y +K*p2z/ERC*K*p2z;
        A21=-K*p2x/ERC*K*p1x -K*p2y/ERC*K*p1y -K*p2z/ERC*K*p1z;
        A12=-K*p1x/ERC*K*p2x -K*p1y/ERC*K*p2y -K*p1z/ERC*K*p2z;
        A11=+K*p1x/ERC*K*p1x +K*p1y/ERC*K*p1y +K*p1z/ERC*K*p1z;
        A=[A11 A12;A21 A22];

        %solve eigenvalue problem
        [~,eigD]=eig(A);
        k02=sort(diag(eigD));
        WN(:,nb)=(crystal.unitcell.a/2/pi)*real(sqrt(k02(1:NSOL)));
        
        %build eigenvalue proble for the empty cell
        A11=+K*p2x/ERMGC*K*p2x +K*p2y/ERMGC*K*p2y +K*p2z/ERMGC*K*p2z;
        A12=-K*p2x/ERMGC*K*p1x -K*p2y/ERMGC*K*p1y -K*p2z/ERMGC*K*p1z;
        A21=-K*p1x/ERMGC*K*p2x -K*p1y/ERMGC*K*p2y -K*p1z/ERMGC*K*p2z;
        A22=+K*p1x/ERMGC*K*p1x +K*p1y/ERMGC*K*p1y +K*p1z/ERMGC*K*p1z;
        A=[A11 A12;A21 A22];
        
        %solve eigenvalue problem for the effective cell
        k02MG=sort(eig(A));
        WNMG(:,nb)=(crystal.unitcell.a/2/pi)*real(sqrt(k02MG(1:NSOL)));

       
        %draw the bands
%         subplot(1,3,[2 3]);
        plot([1:nb],WN(:,1:nb)','Color','b');
        axis tight
        xlim([1 NBETA]);
        ylim([0 wn_max]);
        drawnow;
    end
    
    %% compute the effective refractive index as function of the bloch vector
    absBETA=zeros(NBETA,1);
    for ii=1:NBETA
        absBETA(ii)=norm(BETA(:,ii),2);
    end
    k0fun=(WN(1,:)'*2*pi/crystal.unitcell.a);
    neff=absBETA./k0fun;
    
    %% compute band diagram based on effective refractive index at Gamma symmetry point
    [~,I]=min(absBETA);
    ERLW=ones(size(ER))*neff(I)^2;

    %build convolution matrix
    ERLWC=convmat(ERLW,NP,NQ,NR);

    %main loop - iterate over bloch wave vector
    for nb=1:NBETA

        %calculate wave vector expansion
        Kx=BETA(1,nb)-P*crystal.T1(1)-Q*crystal.T2(1)-R*crystal.T3(1);
        Ky=BETA(2,nb)-P*crystal.T1(2)-Q*crystal.T2(2)-R*crystal.T3(2);
        Kz=BETA(3,nb)-P*crystal.T1(3)-Q*crystal.T2(3)-R*crystal.T3(3);
        K=sqrt(abs(Kx).^2+abs(Ky).^2+abs(Kz).^2);

        %calculate perpendicular polarization vectors
        p1x=speye(NH,NH); p2x=speye(NH,NH);
        p1y=speye(NH,NH); p2y=speye(NH,NH);
        p1z=speye(NH,NH); p2z=speye(NH,NH);
        for m=1:NH
            k=[Kx(m);Ky(m);Kz(m)];
            if norm(k)<1e-10
                p1=[1;0;0];
                p2=[0;1;0];
            else
                v=[4*k(2); 2*k(3); 3*k(1)];
                p1=cross(k,v);
                p1=p1/norm(p1);
                p2=cross(k,p1);
                p2=p2/norm(p2);
            end
            p1x(m,m)=p1(1);
            p1y(m,m)=p1(2);
            p1z(m,m)=p1(3);
            p2x(m,m)=p2(1);
            p2y(m,m)=p2(2);
            p2z(m,m)=p2(3);
        end

        % build the eigenvalue problem
        K=diag(sparse(K(:)));
        A11=+K*p2x/ERLWC*K*p2x +K*p2y/ERLWC*K*p2y +K*p2z/ERLWC*K*p2z;
        A12=-K*p2x/ERLWC*K*p1x -K*p2y/ERLWC*K*p1y -K*p2z/ERLWC*K*p1z;
        A21=-K*p1x/ERLWC*K*p2x -K*p1y/ERLWC*K*p2y -K*p1z/ERLWC*K*p2z;
        A22=+K*p1x/ERLWC*K*p1x +K*p1y/ERLWC*K*p1y +K*p1z/ERLWC*K*p1z;
        A=[A11 A12;A21 A22];

        %solve eigenvalue problem
        k02=sort(eig(A));
        WNLW(:,nb)=(crystal.unitcell.a/2/pi)*real(sqrt(k02(1:NSOL)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DRAW BAND DIAGRAM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cla;
    hold on

    %draw vertical lines at key points of symmetry
    for m=1:length(KT)
        x=KT(m)*[1 1];
        y=[0 wn_max];
        line(x,y,'Color','k','LineStyle','--')
    end

    %draw the bands
    plot([1:NBETA],WN,'-b','LineWidth',1);
    plot([1:NBETA],WNLW,'-r','LineWidth',1);
    plot([1:NBETA],WNMG,'Color',[0.5 0.5 0.5],'LineWidth',1);

    %beta axis
    set(gca,'XTick',KT,'XTickLabel',crystal.KL);
    xlabel('Bloch Wave Vector $\vec{\beta}$','Interpreter','LaTex')
    ylabel('Normalized frequency $\omega=a/\lambda_0$','Interpreter','LaTex');
    YT=[0:0.2:wn_max];
    YL={'0'};
    for m=1:length(YT)
        YL{m}=num2str(YT(m),'%3.1f'); 
    end
    set(gca,'YTick',YT,'YTickLabel',YL);
    %set graphics view
    hold off;
    axis tight;
    ylim([0 wn_max]);
    title('Band Diagram')

    %% plot the extracted effective refractive index
    figure; hold all
    plot(1:NBETA,absBETA,'LineWidth',1)
    plot([1:NBETA],k0fun,'LineWidth',1)
    plot(1:NBETA,neff,'LineWidth',1)
    grid on; box on;
    %beta axis
    set(gca,'XTick',KT,'XTickLabel',crystal.KL);
    xlabel('Bloch Wave Vector $\vec{\beta}$','Interpreter','LaTex')
    legend('$\left|\vec{\beta}\right|$','$k_0$','$n_{eff}$','Interpreter','LaTex')
    axis tight;

    %% plot the difference between the light line and unit cell solutions
    figure;  hold on
    plot(1:NBETA,abs(WNLW(1,:)-WN(1,:)),'LineWidth',1)
    plot(1:NBETA,abs(WNMG(1,:)-WN(1,:)),'LineWidth',1)
    box on; grid on
    set(gca,'XTick',KT,'XTickLabel',crystal.KL);
    xlabel('Bloch Wave Vector $\vec{\beta}$','Interpreter','LaTex')
    legend('PWEM','Maxwell Garnett')
    axis tight;
end