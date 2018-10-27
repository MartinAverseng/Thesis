function [ pol_in_chebBase ] = O1( pol_in_chebBase )

degP = length(pol_in_chebBase) - 1;
eigenvals = [log(2)/2; 1./(2*(1:degP))'];
pol_in_chebBase = eigenvals.*pol_in_chebBase;

end

