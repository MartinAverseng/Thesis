function[P] = SCSD_orderEst(rho,eps)
% Estimation of number of sinc used to approximate unity in [rho, pi-rho]
% with less than eps L2 error

P = - log(eps)/(2*sin(rho));



end