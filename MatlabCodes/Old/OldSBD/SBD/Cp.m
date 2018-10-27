function[C] = Cp(rho)

C = 1./(sqrt(pi)*rho.*abs(besselj(1,rho)));

end