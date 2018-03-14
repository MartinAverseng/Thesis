function [ e ] = calc_err(approx,sol,X,K)

T = acos(X);
T = T(end:-1:1);
Delta = diff(T);
T_phi = T(1:end-1) + Delta/2;
X_phi = cos(T_phi);
sol_h = sol(X_phi);
err_h = approx(:) - sol_h(:);
e = sqrt(err_h'*K*err_h);
% e = sqrt(sum(Delta.*(err_h.^2)));

end

