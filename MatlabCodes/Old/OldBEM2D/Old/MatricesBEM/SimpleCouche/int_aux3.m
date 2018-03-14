function [ I ] = int_aux3(t1,t2,t3,t4,a,b,c,d)

if and(abs(t2-pi)<1e-5,abs(t4-pi)<1e-5)
    E = t2 - t1;
    F = t4 - t3;
    I = (1/12*(a*c*E^2+(-2*a*c*F+(4*c*pi+2*d)*a+2*b*c)*E+3*a*c*F^2+((-8*c*pi-4*d)*a-4*b*c)*F+(6*(c*pi+d))*(a*pi+b)))*(E+F)^2*log(E+F)-(1/12*(a*c*E^2+((4*c*pi+2*d)*a+2*b*c)*E+(6*(c*pi+d))*(a*pi+b)))*E^2*log(E)-(1/4)*F*(F*(a*c*F^2+((-(8/3)*pi*c-(4/3)*d)*a-(4/3)*b*c)*F+(2*(c*pi+d))*(a*pi+b))*log(F)+(1/3*((4*a*c*F^2+((-12*c*pi-6*d)*a-6*b*c)*F+(12*(c*pi+d))*(a*pi+b))*log(2)+a*c*E^2+(-(1/2)*a*c*F+(4*c*pi+2*d)*a+2*b*c)*E+(13/3)*a*c*F^2+((-14*c*pi-7*d)*a-7*b*c)*F+(18*(c*pi+d))*(a*pi+b)))*E);
else
    
    I = (2/3*(((pi^2+(t1+t3)*pi-(3/4)*(t1-t3)^2)*c+(2*(pi-(1/2)*t1+t3))*d)*a+(2*((pi+t1-(1/2)*t3)*c+(3/2)*d))*b))*(pi-(1/2)*t1-(1/2)*t3)^2*log(2*pi-t1-t3)-(2/3*(((pi^2+(t2+t3)*pi-(3/4)*(t2-t3)^2)*c+(2*(pi-(1/2)*t2+t3))*d)*a+(2*((pi+t2-(1/2)*t3)*c+(3/2)*d))*b))*(pi-(1/2)*t2-(1/2)*t3)^2*log(-t2+2*pi-t3)-(2/3)*(pi-(1/2)*t1-(1/2)*t4)^2*(((pi^2+(t4+t1)*pi-(3/4)*(t1-t4)^2)*c+2*d*(pi-(1/2)*t1+t4))*a+(2*((pi+t1-(1/2)*t4)*c+(3/2)*d))*b)*log(-t4+2*pi-t1)+(2/3)*(pi-(1/2)*t2-(1/2)*t4)^2*(((pi^2+(t4+t2)*pi-(3/4)*(t2-t4)^2)*c+(2*(pi-(1/2)*t2+t4))*d)*a+(2*((pi+t2-(1/2)*t4)*c+(3/2)*d))*b)*log(-t4-t2+2*pi)-(1/4*(((t3+t4)*a+2*b)*((t1+t2)*c+2*d)*log(2)+((-(2/3)*pi^2+((5/3)*t1+(5/3)*t2+(5/3)*t3+(5/3)*t4)*pi-(1/2)*t1^2+(-(1/2)*t2+(3/4)*t3+(3/4)*t4)*t1-(1/2)*t2^2+((3/4)*t3+(3/4)*t4)*t2-(1/2)*t3^2-(1/2)*t4^2-(1/2)*t3*t4)*c+(8/3*(pi-(1/4)*t1-(1/4)*t2+(7/8)*t3+(7/8)*t4))*d)*a+(8/3*((pi+(7/8)*t1+(7/8)*t2-(1/4)*t3-(1/4)*t4)*c+(9/4)*d))*b))*(t3-t4)*(t1-t2);
    % cf. Maple
end
end

