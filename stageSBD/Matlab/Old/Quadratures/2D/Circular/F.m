function [H,A0,COS,SIN,s] = F(N)

if N==0
    H = 1;
    
else
    
    n = -N:N;
    G = F(N-1);
    H = -1/2*[0 0 G] + 1i*n.*[0 G 0] +1/2*[G 0 0];
    
end

A0 = H(N+1);
COS = (H(N+2:end) + H(N:-1:1));
SIN = 1i*(H(N+2:end) - H(N:-1:1));

s = '';
if real(A0) ~= 0
    s = sprintf('%.5g ',real(A0));
end
if imag(A0) > 0
    s = [s sprintf('+ %.5g i ',imag(A0))];
elseif imag(A0) < 0
    s = [s sprintf('- %.5g i ',-imag(A0))];
end

for k = 1:length(COS)
    if real(COS(k))>0
         s  = [s sprintf('+ %.5g cos(%dx) ',real(COS(k)),k)]; 
    elseif real(COS(k))<0
         s  = [s sprintf('- %.5g cos(%dx) ',real(-COS(k)),k)];
    end
    if imag(COS(k))>0
         s  = [s sprintf('+ %.5gi cos(%dx) ',imag(COS(k)),k)]; 
    elseif imag(COS(k))<0
         s  = [s sprintf('- %.5gi cos(%dx) ',imag(-COS(k)),k)];
    end
    if real(SIN(k))>0
         s  = [s sprintf('+ %.5g sin(%dx) ',real(SIN(k)),k)]; 
    elseif real(SIN(k))<0
         s  = [s sprintf('- %.5g sin(%dx) ',real(-SIN(k)),k)];
    end
    if imag(SIN(k))>0
         s  = [s sprintf('+ %.5gi sin(%dx) ',imag(SIN(k)),k)]; 
    elseif imag(SIN(k))<0
         s  = [s sprintf('- %.5gi sin(%dx) ',imag(-SIN(k)),k)];
    end
   
end



end

