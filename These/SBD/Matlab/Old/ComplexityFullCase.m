%% Cases 
clear all
close all

const = 0.1;
Ns = [500 1e3 5e3 1e4 5e4 1e5 5e5 1e6];
tol1 = 1e-3;
tol2 = 1e-6;

i=0;
j=0;
for tol = [tol1 tol2]
    j = 0;
    i=i+1;
    for N = Ns
        j=j+1;
        [t1(i,j),t2(i,j),t3(i,j)] = solveLaplaceCase(N,tol,const); 
    end
end
% figure
% loglog(Ns,t1(1,:),'b--');
% hold on
% loglog(Ns,t2(1,:),'b-');
% loglog(Ns,t1(2,:),'r--');
% loglog(Ns,t2(2,:),'r-');
save('t1','t1');
save('t2','t2');
save('t3','t3');