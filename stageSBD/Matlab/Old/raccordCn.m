function [func] = raccordCn(ep,N)

vec = derVec(N,ep);
A = matA(N,ep);
alpha = A\ vec;
func = @(x)(0);
for i = 0:(N-1)
    func = @(x)(func(x) + alpha(i+1)*(x/ep).^(2*i));
end

end

