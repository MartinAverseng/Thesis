function[Beta,err] = SCSDcoeff(rho,P)


A = makeA(rho,pi-rho,P);
b = makeb(rho,pi-rho,P);
Beta = A\b;

err = abs(pi - 2*rho - b'*Beta)^2;

end


function[A] = makeA(a,b,P)

pp = 2*(0:P-1)+1;
x = repmat(pp,P,1);
y = x';


I1 = (b-a)/2 - 1./(4*pp).*(sin(2*pp*b) - sin(2*pp*a));

C = x-y + eps;
% C(logical(eye(P))) = 1;
D = x+y;

I2 = 1./C.*(sin(C*b) - sin(C*a)) ...
    - 1./D.*(sin(D*b) - sin(D*a));
I2 = I2/2;
% I2(logical(eye(P))) = 0;

A = I2;% + diag(I1);

end

function[b] = makeb(a,b,P)
    pp = 2*(0:P-1)+1;
    b = 1./pp.*(cos(pp*a) - cos(pp*b));
    b = b(:);
end