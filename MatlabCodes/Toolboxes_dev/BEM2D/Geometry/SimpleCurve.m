classdef SimpleCurve
    % Parametric curve, can be closed but must not intersect itself
    
    properties
        x; % Function of 1 argument (the parameter t)
        y; %
        I; % interval of the values for the parameter t
        closed@logical;
        boundedSide;
    end
    
    methods
        function[this] = SimpleCurve(xx,yy,II,bS)
            this.x = xx;
            this.y = yy;
            this.I = II;
            a = II(1); b = II(2); Ma = [xx(a);yy(a)]; Mb = [xx(b);yy(b)];
            this.closed = norm(Ma-Mb)<1e-11;
            if this.closed
                assert(logical(exist('bS','var')),'this curve is closed. You must pass the boundedSide (value "left" or "right") in argument');
                assert(ismember(bS,{'left','right'}),['This curve is closed. the value used for argument bS is incorrect. Please use one of the choices : "left" or "right". \n' ...
                     'choose left if the bounded component of the plane lies at the left of the curve, and right otherwise.'])
                this.boundedSide = bS;
            else
                this.boundedSide = 'none';
                
            end
        end
        function[] = plot(this,N)
            if nargin == 1
                N = 100;
            end
            II = this.I;
            t = linspace(II(1),II(2),N);
            xt = this.x(t);
            yt = this.y(t);
            plot(xt,yt);
            axis equal
        end
    end
    
end

