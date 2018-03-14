%% Unitary Test besselJroots
% Brutal test
a = norm(besselzero('J',0,10000)-besselJroots(0,10000));
assert(a < 1e-10);

% Test 1
rho1exact = 2.40482555769577;

startFreq = 0; P = 100; rho = besselJroots(startFreq,P);
Jrho = abs(besselj(0,rho));
assert(norm(Jrho,'inf')<1e-12);
assert(abs(rho(1)-rho1exact)<1e-12);

% Test 2

startFreq = (2000-1/4)*pi; P = 101; rho = besselJroots(startFreq,P);
Jrho = abs(besselj(0,rho));
assert(norm(Jrho,'inf')<1e-12);
leftRootNum = round(rho(1)/pi+1/4);
assert(leftRootNum == 2000);
rhoNums = round(rho/pi+1/4);
assert(ismember(1950,rhoNums));
assert(ismember(2050,rhoNums));

% Test 2 bis

startFreq = (2000-1/4)*pi; P = 100; rho = besselJroots(startFreq,P);
Jrho = abs(besselj(0,rho));
assert(norm(Jrho,'inf')<1e-12);
leftRootNum = round(rho(1)/pi+1/4);
assert(leftRootNum == 2000);
rhoNums = round(rho/pi+1/4);
assert(ismember(1950,rhoNums));
assert(ismember(2049,rhoNums));
assert(~ismember(2050,rhoNums));

% Test 3

startFreq = (100-1/4)*pi; P = 205; rho = besselJroots(startFreq,P);
Jrho = abs(besselj(0,rho));

assert(norm(Jrho,'inf')<1e-12);
rhoNums = round(rho/pi+1/4);
assert(rhoNums(1)==100);
rho = sort(rho);
leftRoot = rho(1);
rightRoot = rho(end);
assert(abs(leftRoot - rho1exact)<1e-12);
assert(abs(rightRoot - (P-1/4)*pi)<1e-3);

rho = true;
disp('success')
