function [ res] = quantitett(j,N,a)

res = 0;
for i = 1:N
    if i ==j
    else
        res = res + abs(cos((i-j)*a)/abs(i-j) - cos((i+j)*a)/(i+j) + (-1)^(i+j)/(i+j));
    end
    
end

end

