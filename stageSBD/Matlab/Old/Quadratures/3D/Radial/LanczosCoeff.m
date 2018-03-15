function [Alpha] = LanczosCoeff(k,P)
% Lanczos sigma-parameters for the quadrature

pp= 2*(0:P-1)+1;
Alpha = 4./(pi*pp).*(sinc(pp/(2*P)*pi)).^k;

Alpha = Alpha(:);

end

