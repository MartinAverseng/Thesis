function [ I ] = int_aux1(t1,t2,t3,t4,a,b,c,d)

% Calcule \int_{t1}^{t2}\int_{t3}^{t4} (a\theta + b)(c\theta' +
% d)\ln[|\theta - \theta'|]d\thetad\þheta'


% les intervalles ne doivent pas se croiser
assert(and(t1<t2,t3<t4));
if or(or(and(t3<t1,t1<t4),and(t3<t2,t2<t4)),and(t1<t3,t3<t2))
    error('intervals must not cross (but they can be equal)')
end
if t1>t3
    I = int_aux1(t3,t4,t1,t2,a,b,c,d);
else
    
    assert(and(t1<=t3,t2<=t4));
    assert(or(t1==t3,t2<=t3));
    
    if t1==t3
        assert(t2==t4);
        E = t2-t1;
        I = (1/4)*E^2*(((E+2*t1)*c+2*d)*((E+2*t1)*a+2*b)*log(E)+((-(7/4)*E^2-6*E*t1-6*t1^2)*c-3*d*(E+2*t1))*a-(3*((E+2*t1)*c+2*d))*b);
        %(cf maple).
    elseif t2==t3
        E = t2 - t1;
        F = t4 - t3;
        I = (1/8)*(E+F)^2*(E^2*a*c+(2*a*c*F+(4*c*t1+(8/3)*d)*a+(4/3)*b*c)*E+a*c*F^2+((4*c*t1+(8/3)*d)*a+(4/3)*b*c)*F+(4*(c*t1+d))*(a*t1+b))*log(E+F)-(1/2)*F^2*(E^2*a*c+(a*c*F+(2*c*t1+d)*a+b*c)*E+(1/4)*a*c*F^2+((c*t1+(2/3)*d)*a+(1/3)*b*c)*F+(c*t1+d)*(a*t1+b))*log(F)-(1/8*(E*(E^2*a*c+((4*c*t1+(8/3)*d)*a+(4/3)*b*c)*E+(4*(c*t1+d))*(a*t1+b))*log(E)+7*F*(E^2*a*c+((9/14)*a*c*F+((18/7)*c*t1+(32/21)*d)*a+(22/21)*b*c)*E+(1/7)*a*c*F^2+(((6/7)*c*t1+(2/3)*d)*a+(4/21)*b*c)*F+(12/7*(c*t1+d))*(a*t1+b))))*E;
        %(cf maple).
    else
        E = t2 - t1;
        F = t3 - t2;
        G = t4 - t3;
        I = (1/8*(((E+F+G+2*t1)^2*c+(8/3*(E+F+G+(3/2)*t1))*d)*a+(4/3*((E+F+G+3*t1)*c+3*d))*b))*(E+F+G)^2*log(E+F+G)-(1/2*(((E+(1/2)*F+(1/2)*G+t1)^2*c+d*(E+(2/3)*F+(2/3)*G+t1))*a+((E+(1/3)*F+(1/3)*G+t1)*c+d)*b))*(F+G)^2*log(F+G)-(1/8)*(E+F)^2*(((E+F+2*t1)^2*c+(8/3*(E+F+(3/2)*t1))*d)*a+(4/3*((E+F+3*t1)*c+3*d))*b)*log(E+F)+(1/2)*F^2*(((E+(1/2)*F+t1)^2*c+d*(E+(2/3)*F+t1))*a+((E+(1/3)*F+t1)*c+d)*b)*log(F)-(7/8*((((3/7)*F^2+((9/7)*E+(3/7)*G+(12/7)*t1)*F+E^2+((9/14)*G+(18/7)*t1)*E+(1/7)*G^2+(6/7)*t1*G+(12/7)*t1^2)*c+(32/21*(E+(7/8)*F+(7/16)*G+(9/8)*t1))*d)*a+(22/21*((E+(4/11)*F+(2/11)*G+(18/11)*t1)*c+(18/11)*d))*b))*E*G;
        %(cf maple).
    end
end





end

