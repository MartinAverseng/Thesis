function [ A ] = matLap(n,N)

A = zeros(2*n+1);
for i = 0:2*n
    for j = max(0,i-2):i
        if j ==i-2
            
            A(i+1,j+1) = (i-2)*(i-N);
            
            
        elseif j==i-1
            A(i+1,j+1) = N-1 - 2*(i-1);
            
            
        elseif j==i
            A(i+1,j+1) = 1;
        end
    end
    
    
end

