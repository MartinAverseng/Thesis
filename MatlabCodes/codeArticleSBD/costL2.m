function[f,g] = costL2(U,A,B)
V = exp(1i*2*pi*U);

f = 1/2*(A*V)'*(A*V);

if nargout > 1 % gradient required
    g = real(-1i*2*pi*diag(V)'*(B*A*V));
end

end