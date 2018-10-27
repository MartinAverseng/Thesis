%!git status

codePath = '/home/martin/Thesis/MatlabCodes';
addpath(genpath(fullfile(pwd,'usefulMatlabFunctions')));
pwd
addDevPack('AbstractMatrix')
addDevPack('SBD')
addDevPack('BEM2D')
cd 'OpenArc'
addpath(genpath(pwd))
disp('Continuer de modifier weightedFEspace, pour faire une correction plus rapide.')