function [ A ] = matrixA( a,b,P )

% Compute frequencies
rho = BesselZeros(P,3*pi/4);

% Compute relevant quantities
C = sqrt(2)./(abs(besselj(1,rho(:))));
D = besselj(0,rho(:)*[a b]);
E = besselj(1,rho(:)*[a b]);

% Assemble A
% - extra-diagonal
A = (-a*(C.*rho.*D(:,1))*(C.*E(:,1))'...
    + a*(C.*E(:,1))*(C.*rho.*D(:,1))'...
    +b*(C.*rho.*D(:,2))*(C.*E(:,2))'...
    -b*(C.*E(:,2))*(C.*rho.*D(:,2))')./(repmat(rho',P,1).^2-repmat(rho,1,P).^2);
% -Diagonal
A(1:(P+1):end) = C.^2.*(...
    b^2/2*(D(:,2).^2+E(:,2).^2) - b./rho.*D(:,2).*E(:,2) -...
    (a^2/2*(D(:,1).^2 + E(:,1).^2) - a./rho.*D(:,1).*E(:,1)));

end

