classdef  logSingK < singularKernel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function[this] = logSingK()
            this.id = 'ln';
            this.k = LogKernel;
        end
        function[Nm] = Nmax(~)
            Nm = Inf;
        end
        function[I] = IN(~,X,A,B,N)
            [a,b,d] = singularKernel.parameters_singInt(X,repmat(A,size(X,1),1),repmat(B,size(X,1),1));
            %I = 1/2*logSingK.Hk(a,b,d,N);
            I = 1/2*logSingK.Ik(a,b,d,N);
        end
    end
    methods (Static)
        function[I0] = I0(a,b,d)
            % Returns \int_{a} \ln(y^2+ d^2)dy
            F0 = @(x)(2*(1/2*x.*log(d.^2 + x.^2)-x + d.*atan(x./d)));
            I0 = F0(b) - F0(a);
        end
        function[I1] = I1(a,b,d)
            % Returns \int_{a}^b \ln(y^2+ d^2)ydy
            F1 = @(x)(1/2*(x.^2 + d.^2).*(log(x.^2 + d.^2)-1));
            I1 = F1(b) - F1(a);
        end
        
        
        %         function[H] = Hk(a,b,d,N)
        %             Hlist = logSingK.allHk(a,b,d,N);
        %             H = Hlist(:,N+1);
        %         end
        function[I] = Ik(a,b,d,N)
            Ilist = logSingK.allIk(a,b,d,N);
            I = Ilist(:,N+1);
        end
        %         function[Hlist] = allHk(a,b,d,N)
        %             % Computes the list of the integrals
        %             % \int_{0}^{b-a} ln(d^2 + (y+a)^2).y^k
        %             % for k = 0..N.
        %             Hlist = zeros(size(a,1),N+1);
        %             Hlist(:,1) = logSingK.I0(a,b,d);
        %             Hlist(:,2) = logSingK.I1(a,b,d) - a.*logSingK.I0(a,b,d);
        %             for k=3:N+1
        %                 Hlist(:,k) = 2/(k)*(gamma(k-1) + alpha(k-1).*Hlist(:,k-1) + beta(k-1).*Hlist(:,k-2));
        %             end
        %             function[a_l] = alpha(l)
        %                 a_l = -l*a;
        %             end
        %             function[b_l] = beta(l)
        %                 b_l = (1-l)/2*(d.^2 + a.^2);
        %             end
        %             function[c_l] = gamma(l)
        %                 c_l = 1/2*(b-a).^(l-1).*F(d.^2 + b.^2) ...
        %                     + (l-1)/2*(...
        %                     (d.^2 + a.^2)/(l-1).*(b-a).^(l-1) ...
        %                     + (2*a/l).*(b-a).^l ...
        %                     + 1/(l+1)*(b-a).^(l+1));
        %             end
        %             function[out] = F(x)
        %                 out = x.*log(x) - x;
        %             end
        %         end
        function[Ilist] = allIk(a,b,d,N)
            Delta = b-a;
            Ilist = zeros(size(a,1),N+1);
            Ilist(:,0+1) = 1./Delta.*logSingK.I0(a,b,d);
            Ilist(:,1+1) = 1./Delta.^2.*(logSingK.I1(a,b,d) - a.*logSingK.I0(a,b,d));
            for k=2:N
                Ilist(:,k+1) = 2/(k+1)*(alpha(k).*Ilist(:,k-1 + 1) + beta(k).*Ilist(:,k-2 + 1)...
                    + gamma(k));
            end
            function[out] = alpha(l)
                out = -a*l./Delta;
            end
            function[out] = beta(l)
                out = -(d.^2 + a.^2)*(l-1)./(2*Delta.^2);
            end
            function[out] = gamma(l)
                out = (l-1)/2*A(l) + a*(l-1)./Delta*A(l-1)...
                    + (d.^2 + a.^2)*(l-1)./(2*Delta.^2)*A(l-2)...
                    + F(d.^2 + b.^2)./(2*Delta.^2);
            end
            function[out] = A(l)
                out = 1/(l+1);
            end
            function[out] = F(x)
                out = x.*log(x) - x;
            end
        end
        %         function[Ilist] = allIk(a,b,d,N)
        %             Ilist = zeros(size(a,1),N+1);
        %             Ilist(:,1) = logSingK.I0(a,b,d);
        %             if N>0
        %                 Ilist(:,2) = logSingK.I1(a,b,d);
        %                 for k=3:N+1
        %                     Ilist(:,k) = beta(k-1).*Ilist(:,k-2) + alpha(k-1);
        %                 end
        %             end
        %             function[A] = A(l)
        %                 A = (b.^(l+1) - a.^(l+1))/(l+1);
        %             end
        %             function[Fy] = F(y)
        %                 Fy = (d.^2 + y.^2).*log(d.^2 + y.^2) - (d.^2 + y.^2);
        %             end
        %             function[alph] = alpha(l)
        %                 alph = (l-1)/(l+1)*(d.^2.*A(l-2) + A(l)) ...
        %                     + 1/(l+1)*(b.^(l-1).*F(b) - a.^(l-1).*F(a));
        %             end
        %             function[bet] = beta(l)
        %                 bet = -(l-1)/(l+1)*d.^2;
        %             end
        %         end
        %         function [ I ] = recIk(a,b,d,k)
        %
        %             % Recursive formula :
        %             %let I{k} = \int_{a}^b ln(y^2 + d^2)*y^k then,
        %             %I_{k} = beta_k.I{k-2} + alpha_k with
        %             % beta_k = -(k-1)/(k+1)
        %             % alpha_k = (k-1)/(k+1)*(d.^2*A(k-2) + A(k))
        %             %       + 1/(k+1)*(b.^(k-1).*F(b) - a.^(k-1).*F(a));
        %             % A and F as defined in subfunctions.
        %             % We compute it iteratively to be fast.
        %             I = 0;
        %             cumbeta = 1;
        %             while (k>1)
        %                 alpha_k = (k-1)/(k+1)*(d.^2.*A(k-2) + A(k)) ...
        %                     + 1/(k+1)*(b.^(k-1).*F(b) - a.^(k-1).*F(a));
        %                 beta_k = -(k-1)/(k+1)*d.^2;
        %                 I = I + cumbeta*alpha_k;
        %                 cumbeta = cumbeta*beta_k;
        %                 k = k-2;
        %             end
        %             if k==0
        %                 I = I + cumbeta*logSingK.I0(a,b,d);
        %             else
        %                 I = I + cumbeta*logSingK.I0(a,b,d);
        %             end
        %
        %
        %             function[A] = A(l)
        %                 A = (b.^(l+1) - a.^(l+1))/(l+1);
        %             end
        %             function[Fy] = F(y)
        %                 Fy = (d.^2 + y.^2).*log(d.^2 + y.^2) - (d.^2 + y.^2);
        %             end
        %
        %
        %         end
    end
end

