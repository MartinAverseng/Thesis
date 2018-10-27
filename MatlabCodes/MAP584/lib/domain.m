function dom = domain(dim, label)
% Extract a domain from a mesh
% - dim : dimension of the domain to be extracted
% - label (optional) : label(s) of elements to extract
% F. Alouges (c) 2017

global mesh

dom.dim = dim;
if nargin==1
    dom.label = -1;
else
    dom.label = label;
end
end

        
        