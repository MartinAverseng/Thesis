function [ func ] = chi()

func = @(x)(tanh(-1./(x.*(x-1))));


end

