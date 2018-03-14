function plotSol(m, ef, Uhddl, option)
% tracée d'une fonction éléments finis sur un maillage
% données = m=maillage, ef=structure EF, 
%           Uhddl = valeur de la fonction EF aux ddl

switch option
    case 1 % trimesh
        % on utilise Uhddl aux sommets pour P1 et P2
        if size(ef.ordre)==0 % EF P0
            error('on ne peut pas tracer d EF P0 avec l option 1')
        end
        trimesh(m.triangles,m.vertices(:,1),m.vertices(:,2),Uhddl(1:size(m.vertices,1)) );
        view(3)
%     case 2 % trisurf
%         % Pr chq triangle K, on trace int_K(Uh)/int_K(1)
%         Uhint = ef.u * Uhddl;
%         Ntri = size(m.triangles,1);
%         Nint = length(integ.w)/Ntri;
%         ValUhK = sum(reshape(Uhint.*integ.w,Nint,Ntri))';
%         Unint = ones(size(Uhint));
%         VolK = sum(reshape(Unint.*integ.w,Nint,Ntri))';
%         ValUhK = ValUhK./VolK;
%         trisurf(m.triangles,m.vertices(:,1),m.vertices(:,2),zeros(size(m.vertices,1),1),ValUhK)
%         axis([min(ef.ddl(:,1))  max(ef.ddl(:,1)) min(ef.ddl(:,2)) max(ef.ddl(:,2))])
%         view(2)
%         colorbar
   otherwise
        error('option non implémentée')

end