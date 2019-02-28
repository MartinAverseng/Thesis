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
m = mshDisk(N,1);
m.vtx(:,1) = 2*m.vtx(:,1);
%[m,~] = mmg(m);
%plot(m)
%axis equal
m.vtx(:,3) = sqrt(1.00000001 - m.vtx(:,1).^2/4 - m.vtx(:,2).^2);
[m,val] = mmg(m);
nrm = sqrt(m.vtx(:,1).^2/4 + m.vtx(:,2).^2 + m.vtx(:,3).^2);
m.vtx(:,1) = m.vtx(:,1)./nrm;
m.vtx(:,2) = m.vtx(:,2)./nrm;
m.vtx(:,3) = m.vtx(:,3)./nrm;
m2 = m;
m2.vtx(:,3) = 0;
plot(m2)
axis equal





