%% Test de l'id�e d'int�grales singuli�res:
% On cherhe � calculer \int_{seg} ln|x - y|/omega(y) dy

%close all
a = -1; %b = -0.99;
for b = [-0.9,-0.99,-0.999,-0.9999]
    x = a:(b-a)/100:a + 10*(b-a);
    omega = @(x)(sqrt(1-x.^2));
    for i = 1:length(x)
        i
        
        % Int�grale approch�e par Matlab
        IMatlab(i) = integral(@(y)(log(abs(x(i) - y))./omega(y)),a,b);
        
        % Int�grale approch�e de Fran�ois :
        l = sqrt((b-a)^2 + (omega(b) - omega(a))^2);
        F = @(x)(x.*log(abs(x)) - x);
        IFrancois(i) = l/(b-a)*(F(b-x(i)) - F(a - x(i)));
    end
    %figure
    plot(((IMatlab-IFrancois)/sqrt(b-a)),'r');
    hold on
end
%fprintf('Matlab : %s \n',num2str(IMatlab));
%fprintf('Fran�ois : %s \n',num2str(IFrancois));