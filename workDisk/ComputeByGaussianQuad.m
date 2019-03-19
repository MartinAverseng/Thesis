function[I] = ComputeByGaussianQuad(sA,tA,sB,tB,sC,tC,theta,K,J)

fun = @(s,t)(K(s,t).*J(s,t));

if cos(theta)<0.03
    I = integral2tri(fun,[sA,tA],[sB,tB],[sC,tC]);
    return
end

% Imatlab = integral2tri(fun,[sA,tA],[sB,tB],[sC,tC]);
% return
funTest = @(s,t)(1./sqrt(t.^2 + s.^2*cos(theta)^2));
% 
% plotFunOnTri(0.1*[sA,tA],0.1*[sB,tB],0.1*[sC,tC],@(s,t)(fun(s,t)./funTest(s,t)))
% 


sAp = sA*cos(theta);
sBp = sB*cos(theta);
sCp = sC*cos(theta);

phiAp = atan2(tA,sAp);
phiBp = atan2(tB,sBp);
phiCp = atan2(tC,sCp);

rhoAp = sqrt(sAp^2 + tA^2);
rhoBp = sqrt(sBp^2 + tB^2);
rhoCp = sqrt(sCp^2 + tC^2);
% 
% figure
% plot([sA sB sC sA],[tA tB tC tA])
% hold on
% plot(0,0,'*');
% text(sA,tA,'A');
% text(sB,tB,'B');
% text(sC,tC,'C');


if isInTri([0,0],[sA,tA],[sB,tB],[sC,tC])
    % ABC = OAB + OAC + OBC
    I = 0;
    [rhos,phis,ws1] = polarQuad(rhoAp,phiAp,rhoBp,phiBp);
    sp = rhos.*cos(phis);
    s = sp/cos(theta);
    t = rhos.*sin(phis);
%     plot(s,t,'o');
    I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws1);
    [rhos,phis,ws2] = polarQuad(rhoBp,phiBp,rhoCp,phiCp);
    sp = rhos.*cos(phis);
    s = sp/cos(theta);
    t = rhos.*sin(phis);
%     plot(s,t,'o');
    I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws2);
    [rhos,phis,ws3] = polarQuad(rhoCp,phiCp,rhoAp,phiAp);
    sp = rhos.*cos(phis);
    s = sp/cos(theta);
    t = rhos.*sin(phis);
    I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws3);
%     plot(s,t,'o');
    %     disp(1);
    %     disp(sum(ws1 + ws2 + ws3) - exactIntRm1Tri([sAp,tA],[sBp,tB],[sCp,tC],[0,0]));
%     disp(I - Imatlab);
else
    % Renumbering
    [phi1,i] = min([phiAp,phiBp,phiCp]);
    [phi3,k] = max([phiAp,phiBp,phiCp]);
    if abs(phi3 - phi1) > pi
        if phiAp < 0
            phiAp = phiAp + 2*pi;
        end
        if phiBp < 0
            phiBp = phiBp + 2*pi;
        end
        if phiCp < 0
            phiCp = phiCp + 2*pi;
        end
        [~,i] = min([phiAp,phiBp,phiCp]);
        [~,k] = max([phiAp,phiBp,phiCp]);
    end
    %     phi2 = setdiff([phiAp,phiBp,phiCp],[phi1,phi3]);
    j = setdiff([1,2,3],[i,k]);
    ABCp = [sAp sBp sCp; tA tB tC];
    ABCp = ABCp(:,[i,j,k]);
    rhoABCp = [rhoAp;rhoBp;rhoCp];
    rhoAp = rhoABCp(i); rhoBp = rhoABCp(j); rhoCp = rhoABCp(k);
    phiABCp = [phiAp;phiBp;phiCp];
    phiAp = phiABCp(i); phiBp = phiABCp(j); phiCp = phiABCp(k);
    sAp = ABCp(1,1); tA = ABCp(2,1);
    sBp = ABCp(1,2); tB = ABCp(2,2);
    sCp = ABCp(1,3); tC = ABCp(2,3);
%     figure
%     plot([sAp sBp sCp sAp],[tA tB tC tA])
%     hold on
%     plot(0,0,'*');
%     text(sAp,tA,'Ap');
%     text(sBp,tB,'Bp');
%     text(sCp,tC,'Cp');
    if isInTri([sBp,tB],[sAp,tA],[sCp,tC],[0,0])
        % ABC = OAC - OBA - OBC
        I = 0;
        [rhos,phis,ws1] = polarQuad(rhoAp,phiAp,rhoCp,phiCp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
%         plot(sp,t,'o');
        I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws1);
        [rhos,phis,ws2] = polarQuad(rhoAp,phiAp,rhoBp,phiBp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
%         plot(sp,t,'o');
        I = I - 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws2);
        [rhos,phis,ws3] = polarQuad(rhoBp,phiBp,rhoCp,phiCp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
        I = I - 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws3);
%         plot(sp,t,'o');
%         disp(sum(ws1 + ws2 + ws3) - exactIntRm1Tri([sAp,tA],[sBp,tB],[sCp,tC],[0,0]));
%         disp(I - Imatlab);
    else
        % ABC = OAB + OBC - OAC
        I = 0;
        [rhos,phis,ws1] = polarQuad(rhoAp,phiAp,rhoCp,phiCp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
%         plot(sp,t,'o');
        I = I - 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws1);
        [rhos,phis,ws2] = polarQuad(rhoAp,phiAp,rhoBp,phiBp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
%         plot(sp,t,'o');
        I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws2);
        [rhos,phis,ws3] = polarQuad(rhoBp,phiBp,rhoCp,phiCp);
        sp = rhos.*cos(phis);
        s = sp/cos(theta);
        t = rhos.*sin(phis);
        I = I + 1/cos(theta)*sum(fun(s,t)./funTest(s,t).*ws3);
%         plot(sp,t,'o');
        %         disp(1);
        %         disp(sum(ws1 + ws2 + ws3) - exactIntRm1Tri([sAp,tA],[sBp,tB],[sCp,tC],[0,0]));
%         disp(I - Imatlab);
    end
end


% if abs(I - Imatlab)>0.1
%     disp('coucou');
% else
%     close all;
% end
% disp(I - Imatlab);
% drawnow


return

% figure
% plot([sA sB sC sA],[tA tB tC tA])
% hold on
% plot(0,0,'*');
% [s1,i] = min([sA,sB,sC]);
% [s3,j] = max([sA,sB,sC]);
% s2 = setdiff([sA,sB,sC],[s1,s3]);
% k = setdiff([1,2,3],[i,j]);
% ts = [tA tB tC];
% t1 = ts(i); t3 = ts(j); t2 = ts(k);

% text(s1,t1,'A1');
% text(s2,t2,'A2');
% text(s3,t3,'A3');

% alpha = (t2 - t1)/(s2 - s1);
% beta = t2 - alpha*s2;

% splot = linspace(s1-(s2 - s1)/2,s2+(s2 - s1)/2,2);
% plot(splot,alpha*splot+beta,'--');

% gamma = (t3 - t1)/(s3 - s1);
% delta = t3 - gamma*s3;

% splot = linspace(s1-(s3 - s1)/2,s3+(s3 - s1)/2,2);
% plot(splot,gamma*splot+delta,'--');

% Tri1 [s1 t1] [s2 t2] [s2 gamma s2 + delta]

% plot([s1 s2 s2 s1],[t1 t2 (gamma*s2 + delta) t1])

% New triangle
sAp = sA/cos(theta);
sBp = sB/cos(theta);
sCp = sC/cos(theta);

rhoA2 = sAp^2 + tA^2;
rhoB2 = sBp^2 + tB^2;
rhoC2 = sCp^2 + tC^2;
rhoIn = isInTri([0 0],[sAp,tA],[sBp,tB],[sCp,tC]);
if rhoIn
    rho2Min = 0;
    phiMin = 0;
    phiMax = 2*pi;
else
    phiA = atan2(tA,sAp);
    phiB = atan2(tB,sBp);
    phiC = atan2(tC,sCp);
    phiMin = min([phiA,phiB,phiC]);
    phiMax = max([phiA,phiB,phiC]);
    if abs(phiMax - phiMin) > pi
        if phiA < 0
            phiA = phiA + 2*pi;
        end
        if phiB < 0
            phiB = phiB + 2*pi;
        end
        if phiC < 0
            phiC = phiC + 2*pi;
        end
        phiMin = min([phiA,phiB,phiC]);
        phiMax = max([phiA,phiB,phiC]);
    end
    rho2Min = 0;
end
rho2Max = max([rhoA2,rhoB2,rhoC2]);

rhoSpace = linspace(sqrt(rho2Min),sqrt(rho2Max),200);
rhoSpace = rhoSpace(2:end);
rhoSpace = rhoSpace - (rhoSpace(2) - rhoSpace(1))/2;
drho = rhoSpace(2) - rhoSpace(1);
phiSpace = linspace(phiMin,phiMax,200);
dphi = phiSpace(2) - phiSpace(1);
[rhos,phis] = meshgrid(rhoSpace,phiSpace);
grid = [rhos(:).*cos(phis(:)), rhos(:).*sin(phis(:))];
grid = grid(isInTri(grid,[sAp,tA],[sBp,tB],[sCp,tC]),:);
sp = grid(:,1);
s = sp*cos(theta);
t = grid(:,2);
% plot(s,t,'o');

I = sum(fun(s,t)./funTest(s,t))*drho*dphi;
% disp(I - Imatlab)

    function [rrhos,pphis,ws] = polarQuad(rhoA,phiA,rhoB,phiB)
        if phiA > phiB
            phiB = phiB + 2*pi;
        end
        if phiB - phiA > pi
            phiAS = phiA;
            rhoAS = rhoA;
            phiA = phiB-2*pi;
            rhoA = rhoB;
            phiB = phiAS;
            rhoB = rhoAS;
        end
        N = 5;
        N2 = 25;
        % We assume here that OAB is oriented trigonometrically
        [phi,wphi] = refQuad(phiA,phiB);
        ssA = rhoA*cos(phiA); ttA = rhoA*sin(phiA);
        ssB = rhoB*cos(phiB); ttB = rhoB*sin(phiB);
        rhoM = (ssA*(ttB - ttA) - ttA*(ssB - ssA))./(cos(phi)*(ttB - ttA) - sin(phi)*(ssB - ssA));
        rrhos = zeros(N2,1);
        pphis = zeros(N2,1);
        ws = zeros(N2,1);
        for id = 1:length(phi)
            [rhoid,wid] = refQuad(0,rhoM(id));
            rrhos(((id-1)*N + 1):(id*N),1) = rhoid;
            pphis(((id-1)*N + 1):(id*N),1) = phi(id);
            ws(((id-1)*N + 1):(id*N),1) = wid*wphi(id);
        end
    end

    function[x,w] = refQuad(alpha,beta)
        a = 1/3*sqrt(5 + 2*sqrt(10/7));
        b = 1/3*sqrt(5 - 2*sqrt(10/7));
        w1 = (322-13*sqrt(70))/1800;
        w2 = (322+13*sqrt(70))/1800;
        x = [0.5*(1-a) ; 0.5*(1-b) ; 0.5 ; 0.5*(1+b) ; 0.5*(1+a)] ;
        w = [w1 ; w2 ;64/225 ; w2 ; w1];
        x = alpha + (beta-alpha)*x;
        w = (beta-alpha)*w;
    end


end