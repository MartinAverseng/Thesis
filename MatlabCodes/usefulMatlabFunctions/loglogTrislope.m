function [  ] = loglogTrislope( X,Y )
% Trace le graphe loglog et ajoute le petit triangle qui indique la pente. 

%1°) Estimation de la pente par régression linéaire
% Les données doivent être positives, donc satisfaire la condition
assert(and(min(X>0)==1,min(Y>0)==1));
% pour pouvoir tracer un loglog

Xlog = log(X);
Ylog = log(Y);

A = polyfit(Xlog,Ylog,1);
slope = A(1);
p = A(2);

% Log(y) approx slope*log(X) + p
% Donc y approx x^slope * exp(p);


loglog(X,Y,'x');
hold on
xl = log(xlim)/log(10);
yl = log(ylim)/log(10);
xcentre = mean(xl);
rx = diff(xl)/6;
ry = diff(yl)/9;
t = 10.^(linspace(xcentre,xcentre+rx,2));
p1 = p + ry;
p2 = p + 2*ry;
loglog(X,X.^slope*exp(p),'k--');
loglog(t,t.^slope*exp(p1),'r');
loglog(t,[t(2) t(2)].^slope*exp(p1),'r');
loglog([t(1) t(1)],[t(1) t(2)].^slope*exp(p1),'r');
s = sprintf('p = %.3f',fix(slope*1000)/1000);
text(t(1),t(2)^slope*exp(p2),s)
grid on;
hold off

end

