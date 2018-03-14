function[f] = costLinf(U,A,dist)

V = exp(1i*2*pi*U);
f = max(abs(A*V).*(max(dist)-dist));

end