classdef AbstractMatrix
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mv_prod
        N1
        N2
        values = [];
    end
    
    methods
        function[this] = AbstractMatrix(mv_p,NN1,NN2)
            if nargin == 0
                mv_p = @(x)(0*x);
                NN1 = 0;
                NN2 = 0;
            elseif nargin ==1
                if isa(mv_p,'AbstractMatrix')
                    this.mv_prod = mv_p.mv_prod;
                    this.N1 = mv_p.N1;
                    this.N2 = mv_p.N2;
                    this.values = mv_p.values;
                    return
                else
                    assert(isa(mv_p,'double'))
                    NN1 = size(mv_p,1);
                    NN2 = size(mv_p,2);
                    this.values = mv_p;
                    mv_p = @(x)(mv_p*x);
                end
            end
            this.mv_prod = mv_p;
            this.N1 = NN1;
            this.N2 = NN2;
        end
        function[out] = size(this,dim)
            if nargin==1
                out = [this.N1,this.N2];
            else
                if dim ==1
                    out = this.N1;
                else
                    out = this.N2;
                end
            end
        end
        function[this_ij] = elem(this,i,j)
                z = zeros(size(this,2),1);
            ei = z; ei(i) = 1;
            ej = z; ej(j) = 1;
            this_ij = ej'*(this*ei);
        end
        function[C] = plus(A,B)
            if and(isa(B,'double'),isscalar(B))
                C = A;
                C.mv_prod = @(x)(A.mv_prod(x) + B);
            elseif and(isa(B,'double'),isscalar(B))
                C = plus(B,A);
            else
                assert(isequal(size(A),size(B)),'Error using  + (Abstract) Matrix dimensions must agree.');
                C = A;
                A = AbstractMatrix(A);
                B = AbstractMatrix(B);                
                C.mv_prod = @(x)(A.mv_prod(x) + B.mv_prod(x));
                try C.values = A.values + B.values;
                catch 
                end
            end
        end
        function[C] = uminus(A)
            C = A;
            if isa(A,'AbstractMatrix')
                C.mv_prod = @(x)(-A.mv_prod(x));
            else
                C = -A;
            end
        end
        function[C] = minus(A,B)
            C = plus(A,-B);
        end
        function[C] = times(A,B)
            if or(isequal(size(A),[1,1]),isequal(size(B),[1,1]))
                C = A*B;
            else
                assert(isequal(size(A),size(B)),'Matrix dimensions do not agree');
                C = full(A).*full(B);
            end
        end
        function[C] = mtimes(A,B)
            assert(or(isa(B,'double'),isa(B,'AbstractMatrix')),sprintf(...
                '2nd arg is of type %s but should be of type ''AbstractMatrix'' or ''double''',...
                class(B)));
            assert(or(isa(A,'double'),isa(A,'AbstractMatrix')),sprintf(...
                '1st arg is of type %s but should be of type ''AbstractMatrix'' or ''double''',...
                class(B)));
            assert(or(or(size(A,2)==size(B,1),isscalar(A)),isscalar(B)));
            if isscalar(B)
                C = A;
                C.mv_prod = @(x)(B*(A*x));
            elseif isscalar(A)                
                C = mtimes(B,A);
            elseif size(B,2) == 1
                % B is a vector, the result must be a vector
                C = A.mv_prod(B);
            else
                A = AbstractMatrix(A);
                B = AbstractMatrix(B);
                C = AbstractMatrix(@(x)(A.mv_prod(B.mv_prod(x))),size(A,1),size(B,2));                
            end
            
        end
        function[A_full] = full(this)
            if this.N2 == 1
                A_full = this.mv_prod(1);
            else
                I = eye(this.N2);
                A_full = zeros(this.N1,this.N2);
                for i = 1:this.N2
                    A_full(:,i) = this.mv_prod(I(:,i));
                end
            end
            
        end
        function[A] = real(this)
            A = AbstractMatrix(@(x)(real(this.mv_prod)),this.N1,this.N2);
        end
        function[x] = mldivide(this,b)
            x = gmres(this.mv_prod,b,[],1e-10,size(this,1));
        end
    end
    
end

