function [zs] = besselZeros(N,tol)
% returns the approximation of the N first bessel zeros at tolerance eps

% Number of zeros from which we try to look in the table rather than
% computing them.
Nthres = 100;

if ~exist('tol','var')
    
    tol = 1e-13;
    
end
if tol < 1e-13
    warning('Unable to compute roots with so much precision, default tol chosen instead (1e-13)');
    tol = 1e-13;
end
if N<=Nthres
    % then we compute it, since it will take less time than charge them
    
    initialGuess = (0:N-1)'*pi + 3*pi/4;
    zs = computeZeros(initialGuess,tol);
    
else
    if exist('BesselZeros/bessZs.mat','file')
        a = load('BesselZeros/bessZs.mat');
        bessZs = a.bessZs;
        K = length(bessZs);
        if N<=K
            zs = bessZs(1:N,1);
        else
            initialGuess = [bessZs; (K:N-1)'*pi + 3*pi/4];
            zs = computeZeros(initialGuess,tol);
            bessZs = zs(:); %#ok (saved)
            if ~isdir('BesselZeros')
                mkdir('BesselZeros');
            end
            save('BesselZeros/bessZs','bessZs');
        end
    else
        initialGuess = (0:N-1)'*pi + 3*pi/4;
        zs = computeZeros(initialGuess,tol);
        bessZs = zs(:); %#ok (saved)
        if ~isdir('BesselZeros')
            mkdir('BesselZeros');
        end
        save('BesselZeros/bessZs','bessZs');
    end
end

zs = zs(:);


end


function[approxZeros] = computeZeros(initialGuess,tol)

approxZeros = initialGuess;
reachedTol = false;

J0 = besselj(0,approxZeros);
indexToUpdate = abs(J0)>tol;

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


end

