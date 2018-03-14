function [ G ] = Helmholtz_perturb( k )
% Computes the Kernel G(|x|) = real(Gk(|x|)) - G0(|x|) = - 1/4 Y_0(k|x|) +
% 1/2*pi ln |x|, which is C^{\infty}(\mathbb{R}^2). 

G = Kernel(@func,@der);

function[y] = func(x)
    y = 0*x;
    I = x>1e-10;
    J = x<=1e-10;
    y(I) = -1/4*(bessely(0,k*x(I))) + 1/(2*pi)*log(x(I));
    y(J) = -1/(2*pi)*log(k)+0.018451073777173;
end

function[y] = der(x)
    y = 0*x;
    I = x>1e-10;
    J = x<=1e-10;
    y(I) = 1/4*k*bessely(1,k*x(I)) + 1./(2*pi*x(I));
    y(J) = 0;
end


end

