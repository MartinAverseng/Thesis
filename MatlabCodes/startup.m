!git add -A 
!git commit -m "Daily Matlab commit"
!git push
!git status

codePath = '/home/martin/Thesis/MatlabCodes';
addpath(genpath('usefulMatlabFunctions'));
addDevPack('AbstractMatrix')
addDevPack('SBD')
addDevPack('BEM2D')
cd 'OpenArc'