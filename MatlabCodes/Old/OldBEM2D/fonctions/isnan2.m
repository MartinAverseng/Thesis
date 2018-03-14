function [ b ] = isnan2( arg )

try
    
    b = isnan(arg);
catch
    b = false;
end


end

