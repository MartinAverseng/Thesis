function [ K ] = rigidite_simpleCoucheP1(X)

N = length(X);
T = acos(X);
T = T(end:-1:1);
Ndof = N+1;

A = matA(T,Ndof);
B = matB(T,Ndof);
C = matC(T,Ndof);
D = matD(T,Ndof);

K = -1/(2*pi)*(A + B + C + D);
end


function[A] = matA(T,Ndof)

A = zeros(Nmailles,Ndof);

for i =1:Ndof
    for j=i:Ndof
        if and(i>1,j>1)
            t1 = T(i-1); t2 = T(i); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            A11_ij = int_aux1(t1,t2,t3,t4,a,b,c,d);
        else
            A11_ij = 0;
        end
        if i> 1
            t1 = T(i-1); t2 = T(i); t3 = T(j); t4 = T(j+1);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
            A12_ij = int_aux1(t1,t2,t3,t4,a,b,c,d);
        else
            A12_ij = 0;
        end
        if j>1
            t1 = T(i); t2 = T(i+1); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = -1/hi; b=1+T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            A21_ij = int_aux1(t1,t2,t3,t4,a,b,c,d);
        else
            A21_ij = 0;
        end
        t1 = T(i); t2 = T(i+1); t3 = T(j); t4 = T(j+1);
        hi = t2 - t1; hj = t4 - t3;
        a = -1/hi; b=1+T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
        A22_ij = int_aux1(t1,t2,t3,t4,a,b,c,d);
        A_ij = A11_ij + A12_ij + A21_ij + A22_ij;
        A(i,j) = A_ij;
    end
end

A = symmetrize(A);

end

function[B] = matB(T,Nmailles)
B = zeros(Nmailles,Nmailles);

for i =1:Nmailles
    for j=i:Nmailles
        if and(i>1,j>1)
            t1 = T(i-1); t2 = T(i); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            B11_ij = int_aux2(t1,t2,t3,t4,a,b,c,d);
        else
            B11_ij = 0;
        end
        if i> 1
            t1 = T(i-1); t2 = T(i); t3 = T(j); t4 = T(j+1);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
            B12_ij = int_aux2(t1,t2,t3,t4,a,b,c,d);
        else
            B12_ij = 0;
        end
        if j>1
            t1 = T(i); t2 = T(i+1); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = -1/hi; b=1+T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            B21_ij = int_aux2(t1,t2,t3,t4,a,b,c,d);
        else
            B21_ij = 0;
        end
        t1 = T(i); t2 = T(i+1); t3 = T(j); t4 = T(j+1);
        hi = t2 - t1; hj = t4 - t3;
        a = -1/hi; b=1+T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
        B22_ij = int_aux2(t1,t2,t3,t4,a,b,c,d);
        
        B_ij = B11_ij + B12_ij + B21_ij + B22_ij;
        B(i,j) = B_ij;
    end
end

B = symmetrize(B);

end

function[C] = matC(T,Nmailles)
C = zeros(Nmailles,Nmailles);

for i =1:Nmailles
    for j=i:Nmailles
        if and(i>1,j>1)
            t1 = T(i-1); t2 = T(i); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            C11_ij = int_aux3(t1,t2,t3,t4,a,b,c,d);
        else
            C11_ij = 0;
        end
        if i> 1
            t1 = T(i-1); t2 = T(i); t3 = T(j); t4 = T(j+1);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
            C12_ij = int_aux3(t1,t2,t3,t4,a,b,c,d);
        else
            C12_ij = 0;
        end
        if j>1
            t1 = T(i); t2 = T(i+1); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = -1/hi; b=1+T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            C21_ij = int_aux3(t1,t2,t3,t4,a,b,c,d);
        else
            C21_ij = 0;
        end
        t1 = T(i); t2 = T(i+1); t3 = T(j); t4 = T(j+1);
        hi = t2 - t1; hj = t4 - t3;
        a = -1/hi; b=1+T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
        C22_ij = int_aux3(t1,t2,t3,t4,a,b,c,d);
        
        C_ij = C11_ij + C12_ij + C21_ij + C22_ij;
        C(i,j) = C_ij;
    end
end

C = symmetrize(C);

end


function[D] = matD(T,Nmailles)
    D = zeros(Nmailles,Nmailles);

for i =1:Nmailles
    for j=i:Nmailles
        if and(i>1,j>1)
            t1 = T(i-1); t2 = T(i); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            D11_ij = int_aux4(t1,t2,t3,t4,a,b,c,d);
        else
            D11_ij = 0;
        end
        if i> 1
            t1 = T(i-1); t2 = T(i); t3 = T(j); t4 = T(j+1);
            hi = t2 - t1; hj = t4 - t3;
            a = 1/hi; b=1-T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
            D12_ij = int_aux4(t1,t2,t3,t4,a,b,c,d);
        else
            D12_ij = 0;
        end
        if j>1
            t1 = T(i); t2 = T(i+1); t3 = T(j-1); t4 = T(j);
            hi = t2 - t1; hj = t4 - t3;
            a = -1/hi; b=1+T(i)/hi; c = 1/hj; d = 1 - T(j)/hj;
            D21_ij = int_aux4(t1,t2,t3,t4,a,b,c,d);
        else
            D21_ij = 0;
        end
        t1 = T(i); t2 = T(i+1); t3 = T(j); t4 = T(j+1);
        hi = t2 - t1; hj = t4 - t3;
        a = -1/hi; b=1+T(i)/hi; c = -1/hj; d = 1 + T(j)/hj;
        D22_ij = int_aux4(t1,t2,t3,t4,a,b,c,d);
        
        D_ij = D11_ij + D12_ij + D21_ij + D22_ij;
        D(i,j) = D_ij;
    end
end

D = symmetrize(D);

end







