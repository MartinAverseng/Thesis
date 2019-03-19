function [tf] = isInTri(X,A,B,C)

tf = (AreaTri(X,A,B) + AreaTri(X,B,C) + AreaTri(X,A,C)  <= AreaTri(A,B,C)+1e-8); 


end

