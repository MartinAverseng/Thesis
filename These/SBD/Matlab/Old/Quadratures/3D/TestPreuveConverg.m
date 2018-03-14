% Test du découpage de l'intégrale 

x1 = 0;
x2 = 0;
x3 = 60;
X = sqrt(x1^2 + x2^2 + x3^2);

N = 40;

%% Terme 1 : 

[xi,w] = Gauss_Legendre1D(N,-1,1);
f = @(u)(exp(1i*x3*u).*besselj(0,sqrt(x1^2+x2^2)*sqrt(1-u.^2)));
g = @(u,thet)(exp(1i*sqrt(1-u.^2)*(x1*cos(thet)+x2*sin(thet)))...
    .*(exp(1i*u*x3)*ones(size(thet))));

T1 = abs(sin(X)/X - 0.5* sum(w.*f(xi)));

for n = 1:200
    [xi,w] = Gauss_Legendre1D(n,-1,1);
    T1(n) = abs(sin(X)/X - 0.5* sum(w.*f(xi)));
    
    thet = (0:(2*n-1))*2*pi/(2*n);
    A = repmat(w,1,2*n)*1/(2*n);
    T2(n) = abs(0.5* sum(w.*f(xi)) - 0.5*sum(sum(A.*g(xi,thet))));
end


plot((log(abs(T1))))
hold on
plot((log(abs(T2))))