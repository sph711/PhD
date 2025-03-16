clear all
close all
clc

%units
degrees =pi/180;

%open figure
figure('Color','w');

%% dashboard
%unit cell
a=1;

%gird parameters
Sx=1*a; Sy=1*a; Sz=a;
NRESLO=10;
NRESHI=10;

%% calculate grids
dx=a/NRESLO;
dy=a/NRESLO;
dz=a/NRESLO;

Nx=ceil(Sx/dx); dx=Sx/Nx;
Ny=ceil(Sy/dy); dy=Sy/Ny;
Nz=ceil(Sz/dz); dz=Sz/Nz;

%axis vectors
xa=[0:Nx-1]*dx;    xa=xa-mean(xa);
ya=[0:Ny-1]*dy;    ya=ya-mean(ya);
za=[0:Nz-1]*dz;    za=za-mean(za);

[X,Y,Z]=meshgrid(ya,xa,za);

%high resolution grid
dx2=a/NRESHI;   Nx2=ceil(Sx/dx2);
dy2=a/NRESHI;   Ny2=ceil(Sy/dy2);
dz2=a/NRESHI;   Nz2=ceil(Sz/dz2);

xa2=linspace(xa(1),xa(end),Nx2);
ya2=linspace(ya(1),ya(end),Ny2);
za2=linspace(za(1),za(end),Nz2);

dx2=xa2(2)-xa2(1);
dy2=ya2(2)-ya2(1);
dz2=za2(2)-za2(1);

[X2,Y2,Z2]=meshgrid(ya2,xa2,za2);

%% calculate list of grating vectors
% here in general one would define a unit cell, do the fourier expansion
% and end up with a thousend planar grading that each has to be varried
% individually. leading to a very slow code. 
% Instead what is done here rather than coming up with a good unit cell,
% expand it and reduce the set, just three planar gratings are used and
% whatever unit cell this provides, one should check if this works. it
% turns out if we take the unit cell and calculate its reciprocal lattice
% vectors of that lattice, and make that the choice for the three planar
% gratings, all 14 of the bravais lattice structures can be found with this

%simple cubic unit cell
K1=2*pi/a *[1;0;0];
K2=2*pi/a *[0;1;0];
K3=2*pi/a *[0;0;1];

KLIST=[K1 K2 K3];
CLIST=[1 1 1];

%% perform the fill fraction sweep
Nxu = ceil(a/dx2); dxu = a/Nxu;
Nyu = ceil(a/dy2); dyu = a/Nyu;
Nzu = ceil(a/dz2); dzu = a/Nzu;

xau = [0:Nxu-1]*dxu; xau=xau-mean(xau);
yau = [0:Nyu-1]*dyu; yau=yau-mean(yau);
zau = [0:Nzu-1]*dzu; zau=zau-mean(zau);

[YU,XU,ZU]=meshgrid(yau,xau,zau);
%build the analog unit cell
NK = length(CLIST);
UC = zeros(Nxu,Nyu,Nzu);
for nk=1:NK
    Kx=KLIST(1,nk);
    Ky=KLIST(2,nk);
    Kz=KLIST(3,nk);
    GA=exp(1i*(Kx*XU+Ky*YU+Kz*ZU));
    UC=UC+CLIST(nk)*GA;
end

%clean up numerical noise
UC=real(UC);
UC=UC-min(UC(:));
UC=UC/max(UC(:));

%generate fill fraction sweep
if 1
    %initialize threshold data
    gth_dat=linspace(0,1,20);
    
    %tabluate fill fraction data
    ff_dat=0*gth_dat;
    for n=1:length(gth_dat)
        %generate binary unit cell
        UCB=UC>gth_dat(n);
        
        %calculate the fill fraction
        ff_dat(n)=sum(UCB(:))/(Nxu*Nyu*Nzu);
        
        %show unit cell
        clf;
        subplot(1,2,1);
        UCS = smooth3(UCB);
        s=isosurface(YU,XU,ZU,UCS,0.5);
        h=patch(s,'FaceColor','g');
        s=isocaps(YU,XU,ZU,UCS,0.5);
        h=patch(s,'FaceColor','g');
        view(10,10)
        axis equal tight;
        xlim([-a/2 a/2]);
        ylim([-a/2 a/2]);
        zlim([-a/2 a/2]);
        camlight; lighting phong;
        title('UNITCELL');
        
        subplot(1,2,2);
        plot(gth_dat(1:n),ff_dat(1:n));
        xlim([gth_dat(1) max(gth_dat)]);
        ylim([0 1]);
        xlabel('Threshold value')
        ylabel('Fill Fraction')
        drawnow
    end
    %save to file
%     save ffdat gth_dat ff_dat;

end
