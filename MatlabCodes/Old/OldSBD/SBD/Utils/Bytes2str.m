function str = Bytes2str(NumBytes)
% BYTES2STR Private function to take integer bytes and convert it to
% scale-appropriate size.

scale = floor(log(NumBytes)/log(1024));
switch scale
    case 0
        str = [sprintf('%.0f',NumBytes) ' b'];
    case 1
        str = [sprintf('%.2f',NumBytes/(1024)) ' ko'];
    case 2
        str = [sprintf('%.2f',NumBytes/(1024^2)) ' Mo'];
    case 3
        str = [sprintf('%.2f',NumBytes/(1024^3)) ' Go'];
    case 4
        str = [sprintf('%.2f',NumBytes/(1024^4)) ' To'];
    case -inf
        % Size occasionally returned as zero (eg some Java objects).
        str = 'Not Available';
    otherwise
       str = 'Over a petabyte!!!';
end
end