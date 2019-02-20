clear all
close all

i = 0;
gams = [];
lambs1 = [];
for l = 1:50
    for m = 0:l
        i = i+1;
        gams(end+1) = gamma((l+m+2)/2)*gamma((l-m+2)/2)/(gamma((l+m+1)/2)*gamma((l-m+1)/2));
        lambs1(end +1) = l*(l+1) - m^2;
    end
end

figure

subplot(2,1,1);

plot(4*gams.^2./lambs1);
ylim([1,2]);

title('Spectre de $S_{\omega}^2(\omega \nabla \cdot (\omega \nabla))^{-1}$','Interpreter','latex')

subplot(2,1,2);

plot(lambs1./(4*gams.^2))
ylim([0.5,1]);

title('Spectre de $S_{\omega}^{-2}\omega \nabla \cdot (\omega \nabla)$','Interpreter','latex')

xlabel('$l^2 + l + m + 1$ (comptage lexicographique des couples $l,m$)','Interpreter','latex')
