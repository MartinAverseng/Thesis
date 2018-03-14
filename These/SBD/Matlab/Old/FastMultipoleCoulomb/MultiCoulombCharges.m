%% Compute the coulombian potential from a distribution of discrete sources

Npoints = 100;
x = randn(Npoints,1);
y = randn(Npoints,1);
z = randn(Npoints,1);
f = abs(randn(Npoints,1));

%% Exact version : 

xcol = x(:);
xrow = xcol';
Xdiff = xcol*ones(1,Npoints) - ones(Npoints,1)*xrow;
ycol = x(:);
yrow = ycol';
Ydiff = ycol*ones(1,Npoints) - ones(Npoints,1)*yrow;
zcol = z(:);
zrow = zcol';
Zdiff = zcol*ones(1,Npoints) - ones(Npoints,1)*zrow;

G = 1./sqrt(Xdiff.^2 + Ydiff.^2 + Zdiff.^2);
G(logical(eye(Npoints))) = 0*ones(Npoints,1);

imagesc(G);
pot = G*f(:);
figure; 
plot(pot);

r = sqrt(x.^2 + y.^2 + z.^2);
Gchap = fft(r);
F = fft(f);

pot2 = ifft(F.*Gchap);
hold on;
plot(pot2);
