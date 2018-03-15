% Test 1: 
rho1exact = 0.893576966279168;
assert(abs(nextY0root(0)-rho1exact) < 1e-12);

% Test 2: 
k = 1:10:1000;
for i = 1:length(k)
    R = nextY0root(k(i));
    assert(abs(bessely(0,R))< 1e-12)
    assert(R >=k(i));
end


disp('success')