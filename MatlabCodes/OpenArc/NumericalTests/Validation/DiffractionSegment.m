%% Test case on the Semi circle domain.

Main;
%curve = semicircle();
%curve = parabola();
% curve = ellipse();
curve = unitSegment();
figure
plot(curve);
theta_inc = -pi/2;

tTot = tic;
N = 400;
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',5,'specialQuadSegs',1:meshAdapt.nseg);
t1 = tic;
Op_assemblingOptions = {'tol',1e-6,'a_factor',8};
k = 20;
Sw = singleLayer(Vh,...
    'Op_opt',Op_assemblingOptions,'correcMethod','constantTerm','k',k);
Swgalerk = Sw.galerkine(Vh,'U');

X =  R2toRfunc(@(Z)(Z(:,1)));
Y =  R2toRfunc(@(Z)(Z(:,2)));
X0 = 0;
Y0 = 10;
% planeWave = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
cilyndWave = R2toRfunc(@(Z)(besselh(0,k*sqrt((Z(:,1) - X0).^2 + (Z(:,2) - Y0).^2))));
% l = Vh.secondMember(-planeWave);
l = Vh.secondMember(-cilyndWave);
[lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,l,[],1e-7,N,Swgalerk.concretePart);
figure
%semilogy(1:length(RESVEC1),RESVEC1,'-o');
%hold on
semilogy(1:length(RESVEC2),RESVEC2,'-o');
drawnow;


%% Compute the radiating solution.

% I want the solution in the square 
% -1 <= x <= 1, -1.5 <= y <= 1.5
a = Sw.Aop.a;
close all;
figure
y1 = -1; y2 = 1;
x1 = -2; x2 = 2;
x1tmp = x1; y1tmp = y1;
x2tmp = x1; y2tmp = y1;
step_x = (x2 - x1)/5;
step_y = (y2 - y1)/5;
while x1tmp < x2
    x2tmp = x1tmp + step_x;
    while y1tmp < y2
        y2tmp = y1tmp + step_y;
        x = linspace(x1tmp,x2tmp,400);
        y = linspace(y1tmp,y2tmp,400);
        [X,Y] = meshgrid(x,y);
        Sw.set_X([X(:) Y(:)],'a',a/5);
        valsDiffr = Sw*lambda2;
%         valsInc = planeWave([X(:) Y(:)]);
        valsInc = cilyndWave([X(:) Y(:)]);
        clear Y;
        valsDiffr = reshape(valsDiffr,size(X,1),size(X,2));
        valsInc = reshape(valsInc,size(X,1),size(X,2));
        clear X;
        amplitude = valsInc+valsDiffr;
        clear valsInc ;
        clear valsDiffr;
        
        res = 20*log(abs(amplitude));
        clear amplitude;
        imagesc(x,y,res);
        hold on
        
        y1tmp = y2tmp;
        
        clear res;
        axis equal
        axis xy
        xlim([x1 x2])
        ylim([y1 y2])
        drawnow
    end
    x1tmp = x2tmp;
    y1tmp = y1;
    y2tmp = y1;
end


colormap gray
caxis([-150 0])
tTot = toc(tTot);

%% Animation.

figure
animateWave(x1,x2,k,amplitude)


