function [ out ] = singKernel(id)

switch id
    case 'ln'
        out = logSingK;
    otherwise
        error('unknown kernel');
end


end

