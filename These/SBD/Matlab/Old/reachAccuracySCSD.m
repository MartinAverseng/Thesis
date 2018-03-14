function[P] = reachAccuracySCSD(rho,Npoint,eps)



    P = 0;
    reached = false;
    
    while and(~reached,P<1000)
        P = P + 1;
        y = QuadVect(SCSDcoeff(rho,P),Npoint,rho,pi-rho);
        err = sqrt(mean((y - 1).^2));
        reached = err < eps; 
    end
    
end
