%% Non-homogenous kernel


close all;
clear all;
clc;

gamma = [1 2 3];
k = 1;
correc = true;
tol = 0;
P = 10;
leg = {};

if correc
    func = @(x)(sin(k*x) + log(x) - 1/4*k*(cos(k) - k*sin(k))*x.^2);
    derivative = @(x)(k*cos(k*x) + 1./x - 1/2*k*(cos(k) - k*sin(k))*x);
else
    func = @(x)(sin(k*x) + log(x));
    derivative = @(x)(k*cos(k*x) + 1./x );
end

kernel = Kernel(func,derivative);

for i = 1:length(gamma)
    gamma_i = gamma(i);
    a_i = gamma_i/P;
    out  = kernel.radialQuadKernel(a_i,tol,'Pmax',P);
    show(out,2);
    hold on
    leg{end+1} = sprintf('$\\gamma = %s$',num2str(gamma_i));
    title('');
    
end

ylim([10^(-9) 1])
set(gca,'XTick',[gamma/P])
ylabel('Error magnitude (log scale)')
xlabel('$r$ (log scale)');

legend(leg,'Location','southwest')
legend boxoff

currentDir = fileparts(mfilename('fullpath'));

if correc
    path = fullfile(currentDir,'ArbitraryKernelGoodSit.tex');
else
    path = fullfile(currentDir,'ArbitraryKernelBadSit.tex');
end
matlab2tikz(path,'width','0.8\plotwidth','parseStrings',false,...
    'extraTikzpictureOptions','trim axis left, trim axis right');


