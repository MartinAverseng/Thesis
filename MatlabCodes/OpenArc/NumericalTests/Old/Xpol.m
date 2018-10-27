function[XpolInChebBase] = Xpol(polInChebBase)

nonConstantCoeffs = polInChebBase(2:end);
constantCoeff = polInChebBase(1);
polInsymmetricChebBase = [1/2*flipud(nonConstantCoeffs); constantCoeff; 1/2*nonConstantCoeffs];
XpolInsymmetricChebBase = 1/2*( [0; 0; polInsymmetricChebBase] + [polInsymmetricChebBase; 0; 0]);
n = (length(XpolInsymmetricChebBase) -1)/2;
coeffNeg = XpolInsymmetricChebBase(n:-1:1);
coeffPos = XpolInsymmetricChebBase(n+2:end);
XpolInChebBase = [XpolInsymmetricChebBase(n+1); coeffNeg + coeffPos ];

end