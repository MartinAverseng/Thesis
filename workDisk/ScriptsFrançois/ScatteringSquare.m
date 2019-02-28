%ScatteringDisk
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clear all
close all

addpath('../Gypsilab/OpenMsh');
addpath('../Gypsilab/OpenMmg');
addpath('../Gypsilab/OpenDom');
addpath('../Gypsilab/OpenFem');
addpath('../Gypsilab/OpenHmx');

N = 10000;
m = mshSquare(N,[1 1]);
x = m.vtx(:,1);
y = m.vtx(:,2);
d = [x+0.5, 0.5-x, y+0.5, 0.5-y];
d = min(d,[],2);
d = 2*sqrt((0.25-x.^2).*(0.25-y.^2));
m.vtx(:,3) = d;
[m,val] = mmg(m);
m2 = m;
m2.vtx(:,3) = 0;
plot(m)
axis equal





