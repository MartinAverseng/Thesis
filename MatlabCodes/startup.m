!git status

codePath = '/home/martin/Thesis/MatlabCodes';
addpath(genpath('usefulMatlabFunctions'));
addDevPack('AbstractMatrix')
addDevPack('SBD')
addDevPack('BEM2D')
cd 'OpenArc'
addpath(genpath(pwd))