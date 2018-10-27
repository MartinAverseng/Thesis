function [ ext ] = polyExtDer( func,a,b,ak,bk,t)

N = (length(ak)-1)/2;
ext1 = zeros(size(t));
ext2 = zeros(size(t));
ext = zeros(size(t));
for i = 1:length(t)
    x = t(i);
    if x <a
        for k = 0:length(ak)-1
            ext1(i) = ext1(i) + ak(k+1)*(x-a).^k/factorial(k);
            ext2(i) = ext2(i) + k*ak(k+1)*(x-a).^(k-1)/factorial(k);
        end
        ext(i) = 2*N*ext1(i).*(x.^(2*N-1)) + ext2(i).*(x.^(2*N));
        
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

