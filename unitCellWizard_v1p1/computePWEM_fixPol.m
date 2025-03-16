function [k02]=computePWEM_fixPol(Q,P,R,T,ERC,BETA,p1)
    
    %number of fourier expansion coefficients
    NH=numel(Q);
    
    %calculate wave vector expansion
    Kx=BETA(1)-P*T(1,1)-Q*T(1,2)-R*T(1,3);
    Ky=BETA(2)-P*T(2,1)-Q*T(2,2)-R*T(2,3);
    Kz=BETA(3)-P*T(3,1)-Q*T(3,2)-R*T(3,3);
    K=sqrt(abs(Kx).^2+abs(Ky).^2+abs(Kz).^2);

    %calculate perpendicular polarization vectors
    p1x=speye(NH,NH); p2x=speye(NH,NH);
    p1y=speye(NH,NH); p2y=speye(NH,NH);
    p1z=speye(NH,NH); p2z=speye(NH,NH);
    for m=1:NH
        k=[Kx(m);Ky(m);Kz(m)];
        p2=cross(k,p1);
        p2=p2/norm(p2)
        
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
end