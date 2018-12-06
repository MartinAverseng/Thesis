%% Test case on some shape with edges.


Main;
addpath(genpath(fullfile('../../DarbasBrunoEtCie/')))
nn = 1000;
curve = unitSegment;
l = length(curve);
k = nn*2*pi/l;
[curve,incWave] = unitSegment(k);

% [curve,incWave,dxf,dyf] = semicircle(k);
%[curve,incWave] = parabola();
%[curve,incWave] = ellipse();
% [curve,incWave] = spirale(k);
% [curve,incWave] = Scurve(k);
% [curve,incWave] = Vcurve(k);
% figure
% plot(curve);


N = fix(20*nn);
meshAdapt = MeshCurve(curve,N,@cos,[-pi,0]);
Vh =  weightedFEspace(meshAdapt,'P1','1/sqrt(1-t^2)',...
    'quadNum',3,'specialQuadSegs',1:meshAdapt.nseg);
M = Vh.Mass.concretePart;
[L,U,P,Q] = lu(M);
invM = @(u)(M\u);%@(u)(Q*(U\(L \(P*u))));
Wh =  weightedFEspace(meshAdapt,'P1','sqrt(1-t^2)',3);
% dM = (Vh.omega_dx_omega)'*AbstractMatrix.spdiag(Vh.W)*(Vh.omega_dx_omega);
% dM = dM.concretePart;
dM =  Wh.dMass.concretePart;
omega2 = Wh.Mass.concretePart;
% K = inv(omega2*dM*omega2);
K1 = dM - k^2*(omega2-M);
K = dM - k^2*omega2;
Np = 40;
theta = pi/3;
keps = k+1i*0.02*k^(1/3);
sqrtDarbasK1 = @(x)(padePrecondDarbas(x,Np,theta,keps,M,K1));
% clear dM;
% clear omega2;
% u = randn(N+1,1);
% plot((dM - k^2*omega2)*u);
% hold on
% plot(real(sqrtDarbasK1(M\sqrtDarbasK1(u))));
Op_opt = {'tol',1e-3,'a_factor',2};
Sw = singleLayer(Vh,...
    'Op_opt',Op_opt,'correcMethod','constantTerm','k',k);
% 
% Nw = hyperSingular_w(Vh,'Op_opt',Op_opt,'k',k);
Swgalerk = Sw.galerkine(Vh,'U');
% T0_scal_phi = Vh.phi'*Vh.W;
% T0_star_galerk = T0_scal_phi*T0_scal_phi'/sum(Vh.W);
% Nwgalerk = Nw.galerkine;
% Prec1 = @(u)(invM(Swgalerk*invM(K*invM(u))));
PrecDarbas = @(u)(invM(sqrtDarbasK1(invM(u))));
% clear M L U Q P dM K K1 omega2 
% Prec2 = @(u)(invM(TrefethenSqrt(dM,4,invM(u),M,1,2*Vh.ndof^2))+ 0*invM(u)+ 1*(1/log(2))^2*invM(T0_star_galerk*invM(u)));
% Prec2 = @(u)(omega2\(TrefethenSqrt(omega2*dMinv*omega2,4,omega2\u,omega2,0.5/Vh.ndof^2,1)));
% Prec3 = @(u)(Nwgalerk*invM(u));
% Prec4 = @(u)(omega2\(Nwgalerk.concretePart*invM(u)));
% Y0 = 0.6;
% X0 = -4;
% cylindWave = R2toRfunc(@(Z)(besselh(0,k*sqrt((Z(:,1) - X0).^2 + (Z(:,2) - Y0).^2))));
% incWave = cylindWave;
secondMemb = Vh.secondMember(-incWave);
% secondMemb = Wh.normalDerivative(dxf,dyf);
% l = Vh.secondMember(-cilyndWave);
% u0 = FE_func(Vh,Swgalerk.concretePart \ l.concretePart);
% t0 = tic;
% [lambda0,FLAG0,RELRES0,ITER0,RESVEC0] = variationalSol(Swgalerk,secondMemb,[],1e-8,N);
% t0 = toc(t0);
% disp(t0)
% t1 = tic;
% [lambda1,FLAG1,RELRES1,ITER1,RESVEC1] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,Prec1);
% t1 = toc(t1);
% disp(t1);
t2 = tic;
[lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,PrecDarbas);
t2 = toc(t2);
disp(t2);
% disp(t2);
% t3 = tic;
% [lambda3,FLAG3,RELRES3,ITER3,RESVEC3] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,Prec3);
% t3 = toc(t3);


% t4 = tic;
% [lambda4,FLAG4,RELRES4,ITER4,RESVEC4] = variationalSol(Swgalerk,l,[],1e-8,N,Prec4);
% t4 = toc(t4);
% % [lambda2,FLAG2,RELRES2,ITER2,RESVEC2] = variationalSol(Swgalerk,l,[],1e-8,N,Swgalerk.concretePart);
% figure
% % semilogy(1:length(RESVEC1),RESVEC1,'-o');
% semilogy(1:length(RESVEC0),RESVEC0/norm(secondMemb.concretePart),'-o');
% hold on;
% % semilogy(1:length(RESVEC1),RESVEC1/norm(Prec1(secondMemb.concretePart)),'-o');
% semilogy(1:length(RESVEC2),RESVEC2/norm(PrecDarbas(secondMemb.concretePart)),'-o');
% semilogy(1:length(RESVEC2),RESVEC2/norm(PrecMartin(secondMemb.concretePart)),'-o');
% semilogy(1:length(RESVEC3),RESVEC3/norm(Prec3(secondMemb.concretePart)),'-o');
% semilogy(1:length(RESVEC4),RESVEC4/norm(Prec4(l.concretePart)),'-o');
% legend({sprintf('Without preconditioner, %s s',num2str(t0)),...
%     sprintf('with preconditioner 1, %s s',num2str(t1)),...
%     sprintf('with preconditioner 2, %s s',num2str(t2)),...
%     sprintf('with preconditioner 3, %s s',num2str(t3)),...
%     sprintf('with preconditioner 4, %s s',num2str(t4))});
% xlabel('Iteration')
% ylabel('Residual')
% set(gca,'FontSize',24);
% drawnow;
    
% [lambda1,~,~,~,RESVEC1,t1] = variationalSol(Swgalerk,secondMemb,[],1e-8,N);
% disp('done')
% [~,~,~,~,RESVEC2,t2] = variationalSol(Swgalerk,secondMemb,[],1e-8,N,Prec);
% figure
% semilogy(1:length(RESVEC1),RESVEC1/norm(secondMemb.concretePart),'-o');
% hold on
% semilogy(1:length(RESVEC2),RESVEC2/norm(Prec(secondMemb.concretePart)),'-o');
% lambda = lambda1;
% clear lambda1
% 
% lambda = variationalSol(Swgalerk,secondMemb,[],1e-8,N,Prec3);
% 
%% Compute the radiating solution.
% a = Sw.Aop.a;
% clear K;
% clear M;
% clear Swgalerk
% clear Wh
% clear L U P Q
% clear secondMemb
% % I want the solution in the square 
% % -1 <= x <= 1, -1.5 <= y <= 1.5
% close all;
% figure
% y1 = -1.5; y2 = 2;
% x1 = -5; x2 = 1.5;
% x1tmp = x1; y1tmp = y1;
% x2tmp = x1; y2tmp = y1;
% step_x = (x2 - x1)/5;
% step_y = (y2 - y1)/5;
% while x1tmp < x2
%     x2tmp = x1tmp + step_x;
%     while y1tmp < y2
%         y2tmp = y1tmp + step_y;
%         x = linspace(x1tmp,x2tmp,600);
%         y = linspace(y1tmp,y2tmp,600);
%         [X,Y] = meshgrid(x,y);
%         Sw.set_X([X(:) Y(:)],'a',a/10);
%         valsDiffr = Sw*lambda2;
% %         valsInc = incWave([X(:) Y(:)]);
%         valsInc = incWave([X(:) Y(:)]);
%         clear Y;
%         valsDiffr = reshape(valsDiffr,size(X,1),size(X,2));
%         valsInc = reshape(valsInc,size(X,1),size(X,2));
%         clear X;
%         amplitude = valsInc+valsDiffr;
%         clear valsInc ;
%         clear valsDiffr;
%         
%         res = 20*log(abs(amplitude))/log(10);
%         clear amplitude;
%         imagesc(x,y,res);
%         hold on
%         
%         y1tmp = y2tmp;
%         
%         clear res;
%         axis equal
%         axis xy
%         xlim([x1 x2])
%         ylim([y1 y2])
%         drawnow
%     end
%     x1tmp = x2tmp;
%     y1tmp = y1;
%     y2tmp = y1;
% end
% 
% % % 
% colormap gray
% caxis([-30 -15]);
% hold on;
% plot(curve);
% axis tight
% 

