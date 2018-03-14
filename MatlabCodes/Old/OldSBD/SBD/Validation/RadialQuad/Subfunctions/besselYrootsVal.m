rho1exact = 0.893576966279168;

% Brutal test
a = norm(besselzero('Y',0,10000)-besselYroots(0,10000));
assert(a < 1e-10);

% Test frequency center functionnality
% Test 1
startFreq = (2000-3/4)*pi; P = 101; rho = besselYroots(startFreq,P);
Yrho = abs(bessely(0,rho));
assert(norm(Yrho,'inf')<1e-12);
leftRootNum = round(rho(1)/pi+3/4);     
assert(leftRootNum == 2000);
rhoNums = round(rho/pi+3/4);
assert(ismember(1950,rhoNums));
assert(ismember(2050,rhoNums));

% Test 2
startFreq = (2000-3/4)*pi; P = 100; rho = besselYroots(startFreq,P);
Yrho = abs(bessely(0,rho));
assert(norm(Yrho,'inf')<1e-12);
leftRootNum = round(rho(1)/pi+3/4);
assert(leftRootNum == 2000);
rhoNums = round(rho/pi+3/4);
assert(ismember(1950,rhoNums));
assert(ismember(2049,rhoNums));
assert(~ismember(2050,rhoNums));

% Test 3
startFreq = (100-3/4)*pi; P = 205; rho = besselYroots(startFreq,P);
Yrho = abs(bessely(0,rho));

assert(norm(Yrho,'inf')<1e-12);
rhoNums = round(rho/pi+3/4);
assert(rhoNums(1)==100);
rho = sort(rho);
leftRoot = rho(1);
rightRoot = rho(end);
assert(abs(leftRoot - rho1exact)<1e-12);
assert(abs(rightRoot - (P-3/4)*pi)<1e-3);

rho = true;
disp('success');