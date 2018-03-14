function [xhat, omega] = formInt(dim, numero)
% [xhat, yhat, omega] = formInt(dim, numero)
% formule d'integration en dimension dim sur Khat
% En 1D Khat = [0, 1]
% En 2D Khat = < A, B, C >
% où A = (0,0) , B = (1,0), C = (0,1)
% numero = numero de la formule d'integration
% xhat : coordonnees des points d'integration
% omega = poids d'intégration correspondants
switch dim
    case 1
        switch numero
            case 1
                % ordre 1 / 1 point
                xhat  = 0.5;
                omega = 1;
            case 2
                % ordre 3 / 2 points
                a = 1/sqrt(3);
                xhat  = [ 0.5 * (1-a) ; 0.5 * (1+a) ];
                omega = [ 0.5         ; 0.5         ];
            case 3
                % ordre 5 / 3 points
                a = sqrt(3/5);
                xhat  = [ 0.5 * (1-a) ; 0.5    ; 0.5 * (1+a) ];
                omega = [ 5/18        ; 8/18 ; 5/18          ];
            case 4
                % ordre 7 / 4 points
                a = sqrt(3/7 - 2/7*sqrt(6/5));
                b = sqrt(3/7 + 2/7*sqrt(6/5));
                c = sqrt(30)/72;
                xhat  = [ (1+a)/2 ; (1-a)/2 ; (1+b)/2 ; (1-b)/2];
                omega = [ 1/4+c   ; 1/4+c   ; 1/4-c   ; 1/4-c  ];
            otherwise
                error('unknown 1D integration formula')
        end
    case 2
        switch numero
            case 1
                % ordre 1 / 1 point
                xhat  = [1/3 1/3];
                omega = 0.5 ;
            case 2
                % ordre 2 / 3 points
                xhat  = [ 0.5 0.5; 0.5 0; 0 0.5];
                omega = [ 1/6 ; 1/6 ; 1/6 ];
            case 3
                % ordre 2 / 3 points
                xhat  = [ 1/6 1/6; 1/6 2/3; 2/3 1/6];
                omega = [ 1/6 ; 1/6 ; 1/6 ];
            case 4
                % ordre 3 / 4 points
                xhat  = [ 1/3 1/3; 1/5 1/5; 3/5 1/5; 1/5 3/5];
                omega = [ -9/32 ; 25/96 ; 25/96 ; 25/96 ];
            case 5
                % ordre 4 / 6 points
                xhat  = [ 1/2 1/2; 1/2 0; 0 1/2; 1/6 1/6; 1/6 2/3; 2/3 1/6];
                omega = [ 1/60 ; 1/60 ; 1/60 ; 3/20 ; 3/20 ; 3/20 ];
            case 6
                % ordre 6 / 7 points
                a = (6 + sqrt(15)) / 21;
                b = (6 - sqrt(15)) / 21;
                A = (155 + sqrt(15)) / 2400;
                B = (155 - sqrt(15)) / 2400;
                xhat  = [ 1/3 1/3; a  a; 1 - 2*a a; a 1-2*a; b b; 1 - 2*b b; b 1-2*b];
                omega = [9/80    ; A   ; A        ; A      ; B  ; B        ; B      ];
            otherwise
                error('unknown 2D integration formula')
        end
end