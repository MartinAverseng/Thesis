function[v,normH10_ab,scal01] = scalProds(a,b,C,kernel,rho)
% This function computes the vector which entries are
% $v(i) = C(i) \int_{a < |x| < b} \nabla f \cdot \nabla g_i$
% with \nabla f given by derivative (assumed to be radial, $|\nabla f| =
% f'(r)$ and g_i is the i_th eigenvector of the Laplace operator with
% Dirichlet b.c., i.e.
% g_i(x) = C(i) J_0(\rho_i |x|)$ with $J_0$ the Bessel function of first kind and 0 order
% C(i) a normalization constant in $H^1_0$ norm and $\rho_i$ the i-th root
% of $J_0$ on the real axis.
% The function uses oscillatory integrals through the matlab function
% "integral" (which is vectorized, but quite slow !).

warning off MATLAB:integral:NonFiniteValue

scalFunc = kernel.scalFunc;
normFunc = kernel.normFunc;

if nargin==0
    run('scalProdsVal.m')
    % Unitary Test
else
    % Vector valued function
    % (i.e. when x is scalar, fun(x) is a vector (here of size
    % [1,length(rho)])
    v = -C(:)'.*scalFunc(a,b,rho);
    
    if a==0
        leftEdge = 0;
    else
        leftEdge = -C(:)'.*scalFunc(0,a,rho);
    end
    
    if b==1
        rightEdge = 0;
    else
        
        rightEdge = -C(:)'.*scalFunc(b,1,rho);
    end
    
    scal01 = (leftEdge + v + rightEdge);
 	normH10_ab = normFunc(a,b);
    
end

v = v(:);
scal01 = scal01(:);

warning on MATLAB:integral:NonFiniteValue

end