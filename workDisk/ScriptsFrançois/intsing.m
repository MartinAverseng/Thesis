%% Test de l'idée d'intégrales singulières:
% On cherhe à calculer \int_{seg} ln|x - y|/omega(y) dy

%close all
a = -1; %b = -0.99;
for b = [-0.9,-0.99,-0.999,-0.9999]
    x = a:(b-a)/100:a + 10*(b-a);
    omega = @(x)(sqrt(1-x.^2));
    for i = 1:length(x)
        i
        
        % Intégrale approchée par Matlab
        IMatlab(i) = integral(@(y)(log(abs(x(i) - y))./omega(y)),a,b);
        
        % Intégrale approchée de François :
        l = sqrt((b-a)^2 + (omega(b) - omega(a))^2);
        F = @(x)(x.*log(abs(x)) - x);
        IFrancois(i) = l/(b-a)*(F(b-x(i)) - F(a - x(i)));
    end
    %figure
    plot(((IMatlab-IFrancois)/sqrt(b-a)),'r');
    hold on
end
%fprintf('Matlab : %s \n',num2str(IMatlab));
%fprintf('François : %s \n',num2str(IFrancois));