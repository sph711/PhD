function C = convmat(A,P,Q,R)
    % CONVMAT           Built convolution matrices for Fourier space methods
    %
    % C = convmat(A,P);         for 1D problems
    % C = convmat(A,P,Q);       for 2D problems
    % C = convmat(A,P,Q,R);     for 3D problems
    %
    % This method constructs a convolution matrix C based on the real space
    % grid structure A

    %% Handle input and outputarguments
    %size of A
    [Nx,Ny,Nz]=size(A);

    %number of spatial harmonics
    if nargin == 2
        Q=1; R=1;
    elseif nargin == 3
        R=1;
    end

    %% Prepare and initialize
    %compute indices of spatial harmonics
    NH=P*Q*R;
    p = [-floor(P/2):+floor(P/2)];
    q = [-floor(Q/2):+floor(Q/2)];
    r = [-floor(R/2):+floor(R/2)];

    %compute fourier coefficients of A
    A = fftshift(fftn(A)) / (Nx*Ny*Nz);

    %compute array indices of center harmonic
    p0 = 1+floor(Nx/2);
    q0 = 1+floor(Ny/2);
    r0 = 1+floor(Nz/2);

    %initialize convolution matrix
    C = zeros(NH,NH);

    %% Calculate convolution matrix
    for rrow = 1:R
    for qrow = 1:Q
    for prow = 1:P
        row = (rrow-1)*Q*P+(qrow-1)*P+prow;
        for rcol = 1 : R
        for qcol = 1 : Q
        for pcol = 1 : P
            col = (rcol-1)*P*Q+(qcol-1)*P+pcol;
            pfft = p(prow) - p(pcol);
            qfft = q(qrow) - q(qcol);
            rfft = r(rrow) - r(rcol);
            C(row,col) = A(p0+pfft,q0+qfft,r0+rfft);
        end
        end
        end
    end
    end
    end
end