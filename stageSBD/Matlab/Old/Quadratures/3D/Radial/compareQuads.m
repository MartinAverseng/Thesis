function[h1,h2] = compareQuads(q,leg,Npoint,rho,eps)


r = linspace(0,pi,Npoint);


h1 = figure;

for k = 1:length(q)
    y{k} = QuadVect(q{k},Npoint,0,pi);
    plot(r,y{k},'DisplayName',leg{k});
    hold on;
end
title('Comparison of quadratures SCSD and Lanczos');
legend show
ylim([1-100*eps,1+50*eps]);
xlim([0,pi]);
set(gca,'XTick',[rho, pi/2,pi - rho]);
set(gca,'XTickLabel',{'\rho', '\pi /2','\pi - \rho'});


h2 = figure;

for k = 1:length(q)
    y{k} = QuadVect(q{k},Npoint,0,pi);
    plot(r,log(abs(y{k}-1)),'DisplayName',leg{k});
    hold on;
end
title('Comparison of quadratures SCSD and Lanczos');
legend show
xlim([0,pi]);
set(gca,'XTick',[rho, pi/2,pi - rho]);
set(gca,'XTickLabel',{'\rho', '\pi /2','\pi - \rho'});



end