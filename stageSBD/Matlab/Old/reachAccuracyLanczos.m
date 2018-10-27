function[P] = reachAccuracyLanczos(rho,kmax,Npoint,eps)

P = zeros(kmax,1);

for k = 1:kmax
    P(k,1) = 0;
    reached = false;
    
    while and(~reached,P(k,1)<1000)
        P(k) = P(k) + 1;
        y = QuadVect(LanczosCoeff(k,P(k)),Npoint,rho,pi-rho);
        err = sqrt(mean((y - 1).^2));
        reached = err < eps; 
    end
    
end
end
