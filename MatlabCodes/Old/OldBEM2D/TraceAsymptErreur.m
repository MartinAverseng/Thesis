%Ce script effectue la validation sur les polynômes de Tchebitchev du
%calcul de l'inverse de S.


Main;

s = 4;
f=@(x)(real(polylog(s,exp(1i*acos(x)))));
sol_theo = @(x)(2*real(polylog(s-1,exp(1i*acos(x)))));
i = 1;
Ns = fix(10.^(1:0.15:2.3));
erreur = 0*Ns;
for j = 1:length(Ns)
    N = Ns(j);
    Xtcheb = cos(linspace(pi,0,N)');
    Delta1 = diff(Xtcheb(:));
    phi_tcheb = Xtcheb(1:end-1) + Delta1/2;
    Xuni = linspace(-1,1,N)';
    T = acos(Xuni);
    phi = T(1:end-1);
    Delta2 = diff(Xuni(:));  
    phi_uni = Xuni(1:N-1) + Delta2/2;
    [lambda,ptsDint,K] = GalerkinP0InvS_tchebNodes(f,N);
    [lambda1,K1] = GalerkinInvS(f,Xtcheb,'P0');
    [lambda2,K2] = GalerkinInvS(f,Xuni,'P0');
    [lambda3,K3] = GalerkinInvS(f,Xtcheb,'P1');
    Delta = pi/N;
    erreur(i) = calc_err(lambda,sol_theo,Xtcheb,K);
    erreur1(i) = calc_err(lambda1,sol_theo,Xtcheb,K1);
    erreur2(i) = calc_err(lambda2,sol_theo,Xuni,K2);
    i = i+1;
end

figure
loglogTrislope(Ns,erreur1);
hold on
loglogTrislope(Ns,erreur2);


title(sprintf('Erreur H^{-1/2}, solution dans H^{%.1f}',s-3/2));
xlabel('Nombre de points du maillage')
ylabel('Norme L^2 de l''erreur')