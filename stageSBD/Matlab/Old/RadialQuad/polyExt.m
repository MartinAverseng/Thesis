function [ ext ] = polyExt( func,a,b,ak,bk,t)

N = (length(ak)-1)/2;
ext = zeros(size(t));
for i = 1:length(t)
    x = t(i);
    if x <a
        for k = 0:length(ak)-1
            ext(i) = ext(i) + ak(k+1)*(x-a).^(k)/factorial(k);
        end
        ext(i) = ext(i).*(x.^(2*N));
        
    elseif x>b
        for k = 0:length(bk)-1
            ext(i) = ext(i) + bk(k+1)*(x-b).^(k)/factorial(k);
        end
        ext(i) = ext(i).*((1-x).^(2*N));
    else
        ext(i) = func(x);
    end
end



end

