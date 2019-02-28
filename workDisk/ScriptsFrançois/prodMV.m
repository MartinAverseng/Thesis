function prodMV(N)

% Definitions
Y = rand(3,N);
V = rand(N,1) + 1i*rand(N,1);
k = 5;

% Initialisations
MV = zeros(N,1) + 1i*zeros(N,1);

% Matrice-Vector product
for i = 1:N
    r     = sqrt( (Y(1,i)-Y(1,:)).^2 + (Y(2,i)-Y(2,:)).^2 + (Y(3,i)-Y(3,:)).^2 );
    MV(i) = ((cos(k*r)+1i*sin(k*r))./(r+0.01)) * V;
end

end