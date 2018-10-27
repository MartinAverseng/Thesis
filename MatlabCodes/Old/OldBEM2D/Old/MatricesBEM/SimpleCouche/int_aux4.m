function [ I ] = int_aux4( t1,t2,t3,t4,a,b,c,d )

phi_i = @(x)(a*x + b);
phi_j = @(y)(c*y + d);
f1 = @(u)(sinc(u/pi)); % attention à la déf de sinc sous matlab...
f2 = @(v)((v<=pi/2).*sinc(v/pi)./iif(v==pi,v,pi-v) + (v>pi/2).*iif(v==0,v*0,sinc((pi-v)/pi)./v));
G = @(x,y)(log(f2((x+y)/2).*f1(abs(x-y)/2)));
f = @(x,y)(phi_i(x).*G(x,y).*phi_j(y));

[xk,wxk] = gauss_lobatto(t1,t2);
[yk,wyk] = gauss_lobatto(t3,t4);

N = length(xk);

X = repmat(xk,1,N);
Y = repmat(yk',N,1);

A = f(X,Y);

I = wxk'*A*wyk;


end

