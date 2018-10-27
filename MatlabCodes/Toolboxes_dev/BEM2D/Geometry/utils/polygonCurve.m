function [ c ] = polygonCurve( vertices,closed,bS)

if ~closed
   if ~exist('bS','var')
       bS = 'none';
   end
end
% Creates the simple curve equal to the polygon with vertices 
if closed
    As = vertices;
    Bs = [vertices(2:end,:) ; vertices(1,:) ];
else
    As = vertices(1:end-1,:);
    Bs = vertices(2:end,:);
end

lengths = sqrt((Bs(:,1) - As(:,1)).^2 + (Bs(:,2) - As(:,2)).^2 );
L = sum(lengths);

%c = SimpleCurve(@(t)(aux(t,1)),@(t)(aux(t,2)),[0 L],bS);
c = SimpleCurve(@(u)(aux(L*(u+1)/2,1)),@(u)(aux(L*(u+1)/2,2)),[-1 1],bS);


    function[d] = aux(t,n)
        M = zeros(length(t),2);
        for j = 1:length(t)
            ij = find(t(j)<=cumsum(lengths),1,'first');
            Lij = lengths(ij);
            A = As(ij,:); B = Bs(ij,:);
            M(j,:) = A + (t(j)-sum(lengths(1:ij-1)))/Lij*(B-A);
        end
        d = M(:,n);
    end

end

