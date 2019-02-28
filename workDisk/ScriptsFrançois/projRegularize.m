function Ms = projRegularize(varargin)

%%% INPUT ANALYSIS
if (length(varargin) == 4)
    X     = varargin{1};
    Ydom  = varargin{2};
    green = varargin{3};
    v     = varargin{4};
else
    Xdom  = varargin{1};
    Ydom  = varargin{2};
    u     = varargin{3};
    green = varargin{4};
    v     = varargin{5};
end

%%% INITIALIZATION
% Mesh data from Y
vtx  = Ydom.msh.vtx;
elt  = Ydom.msh.elt;
ctr  = Ydom.msh.ctr;
nrm  = Ydom.msh.nrm;
stp  = Ydom.msh.stp;
%srf  = Ydom.msh.ndv;

tau  = cell2mat(Ydom.msh.tgt);
nu   = cell2mat(Ydom.msh.nrmEdg);
Nelt = size(elt,1);

%Modification Francois
if ~isempty(Ydom.msh.wgt)
    wgt = Ydom.msh.wgt;
else
    wgt = ones(Nelt,1);
end

% Quadrature data from Y
[Y,Wy,elt2qud] = Ydom.qud;
Yun            = ones(1,size(elt2qud,2));

% Degrees of freedom from Y
[~,elt2dof] = v.dof;
Nbas        = size(elt2dof,2);

% Quadrature data from X
if (length(varargin) == 5)
    [X,Wx] = Xdom.qud;
end
Nx = size(X,1);

% Rangesearch with max(|edge|)_Y
[Ielt,Relt] = rangesearch(X,ctr,1.1*stp(2));                  %%% DEBUG %%%
Mx          = cell(Nelt,1);


%%% RIGHT INTEGRATION WITH REGULARIZATION
for el = 1:Nelt
    % Triangular data for Y
    Sel  = vtx(elt(el,:),:);
    Nel  = nrm(el,:);
    Tel  = reshape(tau(el,:),3,3)';
    NUel = reshape(nu(el,:),3,3)';
    
    % Local size
    edga = Sel(2,:) - Sel(1,:);
    edgb = Sel(3,:) - Sel(1,:);
    edgc = Sel(3,:) - Sel(2,:);
    rMin = 1.1*max([norm(edga),norm(edgb),norm(edgc)]);       %%% DEBUG %%%
    % Quadratures points in interaction
    Iy = elt2qud(el,:);
    Ix = sort(Ielt{el}(Relt{el}<rMin))';

%    Ix = sort(Ielt{el})';
    
%     % Graphical representation                              %%% DEBUG %%%
%     figure(10)
%     plot(Ydom.msh.sub(el))
%     hold on
%     plot3(X(Ix,1),X(Ix,2),X(Ix,3),'*')
%     axis equal
%     hold off
%     pause(1)
    
    % If interactions
    if ~isempty(Ix)
        %%% CORRECTION WITH SEMI-ANALYTIC INTEGRATION
        % Analytical integration
        [Rm1,rRm1,gradRm1,~] = domSemiAnalyticInt3D(X(Ix,:),Sel,Nel,Tel,NUel);
%                 Rm1(:) = 0; rRm1(:) = 0; gradRm1(:) = 0; gradrRm1(:) = 0;    %%% DEBUG %%%
        


        Rm1 = Rm1 * wgt(el);
        rRm1 = rRm1 * wgt(el);
        
        % Vector yg-x
        Xun = ones(length(Ix),1);
        XY1 = Xun * Y(Iy,1)' - X(Ix,1) * Yun;
        XY2 = Xun * Y(Iy,2)' - X(Ix,2) * Yun;
        XY3 = Xun * Y(Iy,3)' - X(Ix,3) * Yun;
                
        % Distance r = |yg-x|
        Rxy              = sqrt(XY1.^2 + XY2.^2 + XY3.^2);
        Rxym1            = 1./Rxy;       
        Rxym1(Rxy<1e-12) = 0;
        
        % Int_el(1/|r|) - Sum_g 1/|yg-x|
        Rm1 = Rm1 - Rxym1 * Wy(Iy);
%                 norm(Rm1)                                   %%% DEBUG %%%

        % Int_el(r/|r|) - Sum_g (yg-x)/|yg-x|
        rRm1(:,1) = rRm1(:,1) - (XY1 .* Rxym1) * Wy(Iy);
        rRm1(:,2) = rRm1(:,2) - (XY2 .* Rxym1) * Wy(Iy);
        rRm1(:,3) = rRm1(:,3) - (XY3 .* Rxym1) * Wy(Iy);
%                 norm(rRm1)                                  %%% DEBUG %%%
        
        % Nullify V
        V = [];
        

        if strcmp(v.typ,'P0')
            % Correction
            if strcmp(green,'[1/r]') && strcmp(v.opr,'[psi]')
                V = Rm1;
            else
                error('domRegularize3D.m : unavailable case')
            end
            
            
        %%% FINITE ELEMENT P1
        elseif strcmp(v.typ,'P1')
            % For each basis function
            for j = 1:Nbas
                % Next dof
                jp1 = mod(j,3) + 1;
                
                % Height from j
                hj = ((Sel(j,:)-Sel(jp1,:)) * NUel(j,:)');
                
                % Scalar product (x-yk).nuj/hj
                tmp = ( (X(Ix,1)-Sel(jp1,1))*NUel(j,1) + ...
                    (X(Ix,2)-Sel(jp1,2))*NUel(j,2) + ...
                    (X(Ix,3)-Sel(jp1,3))*NUel(j,3) ) ./ hj;
                
                % Correction
                if strcmp(green,'[1/r]') && strcmp(v.opr,'[psi]')
                    V(:,j) = Rm1.*tmp + rRm1*NUel(j,:)'/hj;
                    
                elseif strcmp(green,'[1/r]') && strcmp(v.opr,'n*[psi]')
                    Vx        = Rm1.*tmp + rRm1*NUel(j,:)'/hj;
                    V{1}(:,j) = Vx .* Nel(1);
                    V{2}(:,j) = Vx .* Nel(2);
                    V{3}(:,j) = Vx .* Nel(3);
                    
                elseif strcmp(green,'[1/r]') && strcmp(v.opr,'nxgrad[psi]')
                    NxNUj     = cross(Nel,NUel(j,:));
                    V{1}(:,j) = NxNUj(1)/hj .* Rm1;
                    V{2}(:,j) = NxNUj(2)/hj .* Rm1;
                    V{3}(:,j) = NxNUj(3)/hj .* Rm1;
                    
                elseif strcmp(green,'grady[1/r]') && strcmp(v.opr,'n*[psi]')
                    V(:,j) = tmp .* (gradRm1 * Nel');
                    
                elseif strcmp(green(1:end-1),'grady[1/r]') && strcmp(v.opr,'[psi]')
                    ii     = str2double(green(end));
                    V(:,j) = tmp .* gradRm1(:,ii);
                    
                else
                    error('domRegularize3D.m : unavailable case')
                end
            end
        else
            error('domRegularize3D.m : unavailable case')
        end
        
        % Matrix-Vector product
        I = repmat(Ix,1,Nbas);
        J = repmat(elt2dof(el,:),length(Ix),1);
        if iscell(V)
            Mx{el} = [I(:) J(:) V{1}(:) V{2}(:) V{3}(:)];
        else
            Mx{el} = [I(:) J(:) V(:)];
        end
    end
end


%%% LEFT INTEGRATION AND RIGHT REDUCTION
% Left integration matrix
Ndof = size(v.dof,1);
if (length(varargin) == 4)
    Mu = speye(Nx,Nx);
    Mw = speye(Nx,Nx);
    Ms = sparse(Nx,Ndof);
else
    Mu = u.uqm(Xdom);
    Mw = spdiags(Wx,0,length(Wx),length(Wx));
    Ms = sparse(length(u),Ndof);
end

% Regularization matrix
Mx = double(cell2mat(Mx));
if isempty(Mx)
    if iscell(Mu)
        Mx = [1 1 zeros(1,length(Mu))];
    else
        Mx = [1 1 0];
    end
end

% Left integration
if iscell(Mu)
    for i = 1:length(Mu)
        Ms = Ms + Mu{i}' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,2+i),Nx,Ndof);
    end
else
    Ms = Mu' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,3),Nx,Ndof);
end

% Right reduction
[~,Mv] = v.unk;
Ms     = Ms * Mv;
end
