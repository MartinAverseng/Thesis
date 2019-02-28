function [I] = singIntSpecial(A,B,C,X,F1,F2,dFX,N)
% Calcule l'iintégrale sur le triangle ABC de la fonction Y-> 1/norm(F(Y) - F(X))
% où norm est la norme euclidienne, F : R^2 -> R^2 est une certaine
% fonction et dFX est la différentielle de F en X (une matrice 2x2)

% On écrit 1/|F(X) - F(Y)| = 1/|F'(X).(X - Y)| + G(Y)
% où G est plus régulière.

F = @(x,y)([F1(x,y) F2(x,y)]);
X1 = X(1); X2 = X(2);
FX = F(X1,X2);
FX1 = FX(1); FX2 = FX(2);

fsing = @(x,y)(1./normFX_Y(x,y));
fsing2 = @(x,y)(1./normdFXX_Y(x,y));


% % DEBUG
figure
plotFunOnTri(A,B,C,fsing);
title('Original integrand')
figure
plotFunOnTri(A,B,C,@(x,y)(G(x,y)));
title('Integrand once the principal singularity has been removed')

Ap = dFX*A(:);
Bp = dFX*B(:);
Cp = dFX*C(:);
Xp = dFX*X(:);

intfsingExplicit = abs(1./det(dFX))*exactIntRm1Tri(Ap,Bp,Cp,Xp);
intMatlab = integral2tri(fsing2,A,B,C);
freg = @(y1,y2)(G(y1,y2));
intGapprox =  gaussQuadTri(freg,A,B,C,N);
inMatlab =  integral2tri(freg,A,B,C);
I = intfsingExplicit + intGapprox;

    function[out] = G(y1,y2)
        Y = [y1 y2];
        
        if isequal(X,Y)
            out = G(y1+1e-8,y2+1e-8);
        else
            out = fsing(y1,y2) - 1./normdFXX_Y(y1,y2);
        end
    end

    function[out] = normFX_Y(x,y)
        out = sqrt((F1(X1,X2) - F1(x,y)).^2 + (F2(X1,X2) - F2(x,y)).^2);
    end
    function[out] = normdFXX_Y(y1,y2)
        X_Y1 = y1 - X1; X_Y2 = y2 - X2;
        dF1XX_Y = dFX(1,1)*X_Y1  + dFX(1,2)*X_Y2;
        dF2XX_Y = dFX(2,1)*X_Y1  + dFX(2,2)*X_Y2;
        out = sqrt(dF1XX_Y.^2 + dF2XX_Y.^2);
    end
end

