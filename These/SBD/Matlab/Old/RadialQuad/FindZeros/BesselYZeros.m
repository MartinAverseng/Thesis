function [zs] = BesselYZeros(c,k,N,freqCenter)

if ~exist('freqCenter','var')
    freqCenter = 0;
end

xmin = max(freqCenter - pi*N,0);
xmax = freqCenter + 2*pi*N; % Using zn ~ pi*n
if c == Inf
    % Returns the N first zeros of the function
    % Jk(x)
    zs = AllZeros(@(x)(bessely(k,x)),xmin,xmax,20*(N+freqCenter));
    
else
    % Returns the N first zeros of the function
    % c Jk(x) + xJk'(x)
    
    
    switch N
        case 0
            derBess = @(x)(-besselj(1,x));
        otherwise
            derBess = @(x)(0.5*(bessely(k-1,x)-bessely(k+1,x)));
    end
    
    zs = AllZeros(@(x)(c*bessely(k,x)+x*derBess(x)),xmin,xmax,20*(N+freqCenter));
    
    
    
end

[~,I] = sort(abs(zs-freqCenter));
zs = zs(I);

if zs(1)==0
    zs = zs(2:(N+1));
else
    zs = zs(1:N);
end

zs = zs(:);

end