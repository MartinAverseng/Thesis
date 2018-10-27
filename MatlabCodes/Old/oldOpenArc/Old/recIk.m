function [ I ] = recIk(a,b,d,k)

% Recursive formula : 
%let I{k} = \int_{a}^b ln(y^2 + d^2)*y^k then, 
%I_{k} = beta_k.I{k-2} + alpha_k with
% beta_k = -(k-1)/(k+1)
% alpha_k = (k-1)/(k+1)*(d.^2*A(k-2) + A(k)) 
%       + 1/(k+1)*(b.^(k-1).*F(b) - a.^(k-1).*F(a));
% A and F as defined in subfunctions. 
% We compute it iteratively to be fast. 
I = 0;
cumbeta = 1;
while (k>1)
    alpha_k = (k-1)/(k+1)*(d.^2.*A(k-2) + A(k)) ...
        + 1/(k+1)*(b.^(k-1).*F(b) - a.^(k-1).*F(a));
    beta_k = -(k-1)/(k+1)*d.^2;
    I = I + cumbeta*alpha_k;
    cumbeta = cumbeta*beta_k;
    k = k-2;
end
if k==0
    F0 = @(x)(2*(1/2*x.*log(d.^2 + x.^2)-x + d.*atan(x./d)));
    I0 = F0(b) - F0(a);
    I = I + cumbeta*I0;
else
    assert(k==1);
    F1 = @(x)(1/2*(x.^2 + d.^2).*(log(x.^2 + d.^2)-1));
    I1 = F1(b) - F1(a);
    I = I + cumbeta*I1;
end


    function[A] = A(l)
        A = (b.^(l+1) - a.^(l+1))/(l+1);
    end
    function[Fy] = F(y)
        Fy = (d.^2 + y.^2).*log(d.^2 + y.^2) - (d.^2 + y.^2);
    end


end

