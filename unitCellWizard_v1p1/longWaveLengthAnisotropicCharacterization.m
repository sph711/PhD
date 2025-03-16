function longWaveLengthAnisotropicCharacterization(crystal,pwemParam,sh,th)
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
    FF=1-sum(sum(sum(UC)))/prod(size(UC));
    
    %compute effective unit cell permittivity according to maxwell garnett
    epEffMG1=pwemParam.nH^2*(2*FF*(pwemParam.nL^2-pwemParam.nH^2)+pwemParam.nL^2+2*pwemParam.nH^2)/(2*pwemParam.nH^2+pwemParam.nL^2-FF*(pwemParam.nL^2-pwemParam.nH^2));
    %convert unit cell to real materials
    ER=pwemParam.nL^2 +(pwemParam.nH^2-pwemParam.nL^2)*UC;
    
    %% compute plane wave expansion parameters
    NP=pwemParam.NP; NQ=pwemParam.NQ; NR=pwemParam.NR;
    
    %build convolution matrix
    ERC=convmat(ER,NP,NQ,NR);
    
    %number of spatial harmonics
    NH=NQ*NP*NR;

    %calculate fourier expansion coefficiens
    p=[-floor(NP/2):floor(NP/2)];
    q=[-floor(NQ/2):floor(NQ/2)];
    r=[-floor(NR/2):floor(NR/2)];
    [Q,P,R]=meshgrid(q,p,r);
    
    %% compute eigensolutions 
    betaA=0.5/1;
    ep=zeros(3,2);
    %propagation x 
    beta=betaA*[1;0;0]; 
    k02=computePWEM(Q,P,R,[crystal.T1,crystal.T2,crystal.T3],ERC,beta);
    ep(1,1)=(betaA/real(sqrt(k02(1)))).^2;
    ep(1,2)=(betaA/real(sqrt(k02(2)))).^2;
    %propagation y 
    beta=betaA*[0;1;0];
    k02=computePWEM(Q,P,R,[crystal.T1,crystal.T2,crystal.T3],ERC,beta);
    ep(2,1)=(betaA/real(sqrt(k02(1)))).^2;
    ep(2,2)=(betaA/real(sqrt(k02(2)))).^2;
    %propagation z 
    beta=betaA*[0;0;1];
    k02=computePWEM(Q,P,R,[crystal.T1,crystal.T2,crystal.T3],ERC,beta);
    ep(3,1)=(betaA/real(sqrt(k02(1)))).^2;
    ep(3,2)=(betaA/real(sqrt(k02(2)))).^2;
    ep
    
end