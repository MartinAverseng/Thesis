function [res] = Oi(i,pol_in_chebBase)
    if i == 1
        res = O1(pol_in_chebBase);
    else
        res = Xpol(Oi(i-1,pol_in_chebBase)) - Oi(i-1,Xpol(pol_in_chebBase));
    end
end