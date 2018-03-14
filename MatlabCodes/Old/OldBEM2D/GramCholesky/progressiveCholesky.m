function [P, alpha] = progressiveCholesky( Afunc,bfunc,tolmonitor,tol,batch_size)
% Input an infinite matrix A, and an infinite vector p. For each value of
% P, denote by A_P='afunc(1,P)' and b_P = 'bfunc(1,P)' the extracted PxP and P 
% first terms of the latter and alpha_P the solution of the linear system 
% A_P alpha_P = b_P. The tolerance monitor 'tolmonitor' maps a vector
% alpha_P with an error 'tolmonitor(alpha_P)'. 
% While the error is not below the required 'tol', a higher value of P is 
% chosen and the solution alpha_P updated. At each step, P is increased by 
% batch_size. 
% The matrix A_P is update as A_P = [A_P,zeros(batch_size);zeros(...);zeros(...)] +
% Afunc(P+1,P + batch_size);

%% Initialization
alpha = [];
beta = [];
P = 0;
b = [];
A = [];
U = [];

%% Main loop
while tolmonitor(alpha) > tol
    A11 = A;
    A = [A, zeros(batch_size);zeror(batch_size),zeros(batch_size)] + Afunc(P+1,P + batch_size);
    A22 = A(P+1:P+batch_size);
    A21= A(P+1:P+batch_size,1:P);
    U1 = U;
    U2 = chol(A2);
    P = P+batch_size;
    beta1 = beta;
    b2 = bfunc(P+1,P+batch_size);
    beta2 = U2'\(b2 - A21*alpha);
    alpha2 = U2\beta2;
    alpha1 = alpha - U1\(U1'\(A21'*alpha2));
    alpha = [alpha1;alpha2];
    beta = [beta1,beta2];
    U = []
end

