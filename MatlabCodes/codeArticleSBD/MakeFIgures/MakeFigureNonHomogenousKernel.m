%% Non-homogenous kernel


close all;
clear all;
clc;

gamma = linspace(0.5,8,10);
k = 50;

%Delta_n_at_1 = [k*(cos(k) - k*sin(k));...
%    k^4*sin(k)-2*k^3*cos(k)+k^2*sin(k)+k*cos(k);...
%   3*k^5*cos(k)-k^6*sin(k)-6*k^3*cos(k)-3*k^4*sin(k)+9*k^2*sin(k)+9*k*cos(k)];

Delta_n_at_1 = [-2424.568848; 6.084459932*10^6; -1.525010569*10^10];

markers = {'-o','-x','-s','-d'};
tol = 0;
P = 100;
lengendd  = {};
err = zeros(length(gamma),length(Delta_n_at_1)+1);
Linf_vs_gamma = figure;
monitor_fig = figure;
leg = {};
for p = 0:length(Delta_n_at_1)
    close(monitor_fig);
    monitor_fig = figure;
    if p==0
        func = @(x)(cos(k*x)./x);
        derivative = @(x)(-k*sin(k*x)./x - cos(k*x)./x.^2);
        kernel = Kernel(func,derivative);
    else
        Q = multiDirBesselSum(Delta_n_at_1(1:p),besselJroots(k,p)+pi/2);
        dQ = diff(Q);
        func = @(x)(cos(k*x)./x - Q(x));
        derivative = @(x)(-k*sin(k*x)./x - cos(k*x)./x.^2 - dQ(x));       
        kernel = Kernel(func,derivative);
    end
    for i = 1:length(gamma)
        gamma_i = gamma(i);
        a_i = gamma_i/P;
        out  = kernel.radialQuadKernel(a_i,tol,'Pmax',P);
        figure(monitor_fig)
        show(out,2);
        err(i,p+1) = out.errLinf;
        disp(err(i,p+1));
        hold on
        lengendd{end+1} = sprintf('$\\gamma = %s$',num2str(gamma_i));
        title('');
        if p==0
            ylabel('Error magnitude (log scale)')
        else
            set(gca,'YTick',[]);
        end
        set(gca,'XTick',[]);
    end
    
    legend(lengendd);
    xlabel('$r$ (log scale)');
    legend show;
    legend boxoff;
    box off;
    
    figure(Linf_vs_gamma)
    semilogy(gamma,err(:,p+1),markers{p+1});
    leg{end+1} = sprintf('$n = %s$',num2str(p));
    hold on
    
    
end
close(monitor_fig);
figure(Linf_vs_gamma);
legend(leg,'Location','southwest')
legend boxoff
xlabel('$\gamma$');
ylabel('$L^{\infty} error$');
currentDir = fileparts(mfilename('fullpath'));
path = fullfile(currentDir,'LinfVsGamma_arbitraryKern2.tex');
matlab2tikz(path,'width','0.8\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');


