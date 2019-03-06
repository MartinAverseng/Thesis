classdef weightedDom < dom
    
    
    
    properties
        quadPoints = [];
        weights = [];
        elt2qud = [];
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function domain = weightedDom(msh,Ngss)
            % Construction de la quadrature:
            domain.msh = msh;
            j = 0;
            n = size(quadSphTri([0,0],[0,1],[1,0],Ngss),1);
            for i = 1:size(msh.elt,1)
                tri = msh.elt(i,:);
                triCoords = msh.vtx(tri,:);
                A = triCoords(1,:);
                B = triCoords(2,:);
                C = triCoords(3,:);
                [X,W] = quadSphTri(A,B,C,Ngss);
                X = [X zeros(n,1)]; % Recast as 3D points
                domain.quadPoints = [domain.quadPoints; X];
                domain.weights = [domain.weights; W];
                domain.elt2qud(i,:) = [j+1, j+2, j+3];
                j = j+3;
            end
            domain.gss = n;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot(varargin)
            domain = varargin{1};
            spc    = '.b';
            if (nargin == 2)
                spc = varargin{2};
            end
            X = domain.qud;
            plot3(X(:,1),X(:,2),X(:,3),spc)
        end
        
        function plotNrm(varargin)
            domain = varargin{1};
            spc    = 'b';
            if (nargin == 2)
                spc = varargin{2};
            end
            Xqud = domain.qud;
            Vnrm = domain.qudNrm;
            quiver3(Xqud(:,1),Xqud(:,2),Xqud(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),spc);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LENGTH
        function l = length(domain)
            l = size(domain.qud,1);
        end
        
        % SIZE
        function s = size(varargin)
            s = size(varargin{1}.qud);
            if (nargin == 2)
                s = s(varargin{2});
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GAUSSIAN QUADRATURE
        function [X,W,elt2qud] = qud(domain)
            X = domain.quadPoints;
            W = domain.weights;
            elt2qud = domain.elt2qud;
        end
        
        % GAUSSIAN NORMALES
        function N = qudNrm(domain)
            x = domReference(domain);
            N = zeros(size(x,1)*size(domain.msh.elt,1),size(domain.msh.vtx,2));
            for j = 1:size(x,1)
                idx      = (j:size(x,1):size(x,1) * size(domain.msh.elt,1))';
                N(idx,:) = domain.msh.nrm;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INTEGRAL FCT
        function I = integral(varargin)
            switch nargin
                case 2
                    I = domIntegral2(varargin);
                case 3
                    I = domIntegral3(varargin);
                case 4
                    I = domIntegral4(varargin);
                case 5
                    I = domIntegral5(varargin);
                case 6
                    I = domIntegral6(varargin);
                case 7
                    I = domIntegral7(varargin);
                case 8
                    I = domIntegral8(varargin);
                otherwise
                    error('dom.m : unavailable case')
            end
        end
        
        % SINGULAR REGULARIZAION
        function S = regularize(varargin)
            mesh = varargin{end}.msh;
            if (size(mesh,2) == 3)
                S = domRegularize3D(varargin);
            elseif (size(mesh,2) == 2) && is2d(mesh)
                S = domRegularize2D(varargin);
            end
        end
        
        % INTERPOLATION
        function I = interpolate(varargin)
            domain = varargin{1};
            u      = varargin{2};
            v      = varargin{3};
            M      = integral(domain,u,u);
            Fv     = integral(domain,u,v);
            if (nargin == 4)
                f = varargin{4};
            else
                f = 1;
            end
            if iscell(Fv)
                I{1} = M \ (Fv{1} * f);
                I{2} = M \ (Fv{2} * f);
                I{3} = M \ (Fv{3} * f);
            else
                I = M \ (Fv * f);
            end
        end
        
        % L2 AND H1 ERRORS
        function err = diff(domain, fe, sol, ref, type )
            err = domDifference(domain, fe, sol, ref, type );
        end
    end
end
