function [ phi ] = padePrecondDarbas( l,Np,theta,keps,I,D)
% g is the Dirichlet data on the boundary, 

k = real(keps);
[ C0,Aj,Bj ] = rotatingPadeRacine(Np,theta);
phi = I\(1i*k*C0*l);

for j = 1:Np
    phi_j = (-Bj(j)/keps^2*D + I)\l;
    phi = phi + (-1i*k*Aj(j)/keps^2)*(I\(D*phi_j));
end


end

