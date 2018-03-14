function [ W ] = quadWeights( X )

ns = 0:length(X)-1;
aux = ((1+(-1).^(ns))./(1+ns))';
VanderMonde = fliplr(vander(X));
W = VanderMonde\aux;

end

