%% Compute and store bessel function zeros :

N = 50; % N first zeros 
tol = 1e-13; 
% I don't manage to be more precise...

approxZeros = (0:N-1)*pi + 3*pi/4;
reachedTol = false;

J0 = besselj(0,approxZeros);
indexToUpdate = 1:N;

while (~reachedTol)
    appZ = approxZeros(indexToUpdate);
    J0 = J0(indexToUpdate);
    J1 = besselj(1,appZ);
    
    
    appZNew = appZ + J0./J1;
    approxZeros(indexToUpdate) = appZNew;
    
    J0 = besselj(0,approxZeros);
    indexToUpdate = abs(J0)>tol;
    reachedTol = sum(indexToUpdate)==0; 
    
end

bessZs = approxZeros;
save('BesselZeros/bessZs','bessZs');

