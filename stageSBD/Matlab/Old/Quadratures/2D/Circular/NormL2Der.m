%% Norme L2 de la dérivée n-ième

clear all;

Ns = 1:3000;
H = 1;
vec = [-Ns 0 Ns];

A = diag([-1i*-fliplr(Ns) 0 -1i*Ns]);
    
for i = 1:length(vec)
    if i - 1 >= 1
        A(i-1,i) = -1/2;
    end
    if i+1 <= size(A,1)
        A(i+1,i) = 1/2;
    end
end

x0 = 1;
x = x0;
trueX = x0;
fact = 0;

for n = Ns
    nvec1 = (1:(max(Ns)-n))';
    vec = max(Ns)+1+(-n:n)';
    nvec2 = ((max(Ns)+1+n+1):(2*max(Ns)+1))';
    fact = fact + log(max(abs(x)));
    x = A(vec,vec)*[0;x;0]/abs(max(x));
    logNorm(n) = 0.5*log(sum(abs(x).^2)) + fact;
    if n<=15
        trueX = A(vec,vec)^n*[zeros((length(vec)-1)/2,1); ...
            x0; ...
            zeros((length(vec)-1)/2,1)];
        trueNorm(n) = sqrt(sum(abs(trueX).^2));
        assert(abs(trueNorm(n) - exp(logNorm(n)))<=1e-6)
        assert(norm(trueX - x*exp(fact)) <= 1e-6);
    end
    
end

