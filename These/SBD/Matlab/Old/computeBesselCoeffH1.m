function [alpha,rho,w0,resH1square] = computeBesselCoeffH1(a,b,func,derivative,tolH1,tolInf,askGraph)


if ~exist('askGraph','var')
    askGraph = false;
end


%% Defining parameters, initialize
timeTotal = tic;
batch_size = min(max(fix((-log(tolH1)+1)/a/3)+1,100),500);

reachedTol = false;
rho = [];
nextFreqApprox = 3*pi/4; % Approximate of the next frequency to be computed
P = 0; % Number of components to compute before the start of the next for loop
A = []; % Matrix of general term (ei|ej)_{L^2[ring]}
B = []; % Matrix for which newE_i = sum(Bij oldE_j)
BBprime = []; % Matrice stoquée pour ne pas recalculer les B*B' en entier
normSQ = [];
beta = [];
alpha = [];

totalOrtho = 0; % Time spent in orthonormalization
totalBessel = 0; % Time spent evaluating bessel functions
totalLinfNorm = 0; % Time spent in evaluation of Linf norm
    
%% Enter in the main loop
while ~reachedTol
    % Pre-allocate rho :
    rho_old = rho;
    rho = zeros(P+batch_size,1);
    rho(1:P) = rho_old;
    

    % Compute next values of rho
    tFreq = tic;
    f  = @(r) besselj(0,r);
    df = @(r) - besselj(1,r);
    init_guess = nextFreqApprox + (0:batch_size-1)*pi;
    rho(P+1:end) = newton(f,df,init_guess);
    tFreq = toc(tFreq);
    fprintf(sprintf('%d frequencies computed in %5f seconds \n',batch_size,tFreq));
    %     % Check if only zeros
    %
    %     if max(abs(f(rho)))>1e-12
    %         error('Problem in newton search for zeros')
    %     end
    %     % Check if all zeros
    %     figure
    %     t = linspace(0,rho(end),15*length(rho));
    %     plot(t,besselj(0,t));
    %     hold on
    %     plot(rho,besselj(0,rho),'*r')
    %     % Check complete
    %
    
        %% Values to test for exiting the loop with Linf criterion
    t1 = tic;
    tTestLinf = linspace(a,2*a,fix(100*batch_size*a)+1);  % 50 points per freq
    C = sqrt(2)./(abs(besselj(1,rho(:))));
    D = besselj(0,rho(:)*[a b]);
    E = besselj(1,rho(:)*[a b]);
    storedValsJ0Mat = diag((C./rho))*besselj(0,rho(:)*tTestLinf(:)');
    storedValsFunc = func(tTestLinf);
    
    
    
    %% Compute the next scalar products of func with basis functions
    % Quadrature to compute the scalar product of func against bessel functions
    resolNGL = 5;
    NGL = min(resolNGL*length(rho),4000);
    [x,w] = Gauss_Legendre1D(NGL,a,b);
    ff = w.*x.*derivative(x);
    
    JJ = diag(C)*-besselj(1,rho(:)*x(:)');
    t1 = toc(t1);
    totalBessel = totalBessel + t1;
    fprintf('Evaluation of J0 and J1 in %d seconds \n',t1);
    f_ei = JJ*ff;
    resH1square = sum(w.*x.*derivative(x).^2)-sum(beta.^2);
    
    
    %% pre/re-allocate all vectors and matrix
    Aold = A; Bold = B; BBprime_old  = BBprime; normSQold = normSQ; betaOld = beta; alphaOld = alpha;
    A = zeros(P+batch_size); B = zeros(P+batch_size); BBprime = zeros(P+batch_size); normSQ = zeros(P+batch_size,1); beta = zeros(P+batch_size,1); alpha = zeros(P+batch_size,1);
    A(1:P,1:P) = Aold; B(1:P,1:P) = Bold; BBprime(1:P,1:P) = BBprime_old;  normSQ(1:P,1) = normSQold; beta(1:P,1) = betaOld; alpha(1:P,1) = alphaOld;
    
    %% Start Gram-Schmidt

    tOrtho = tic;
    for i = P+1:P+batch_size
        % Computing A
        for j = 1:i
            if i==j
                J0a = D(i,1);
                J1a = E(i,1);
                J0b = D(i,2);
                J1b = E(i,2);
                A(i,j) = C(i)^2*(...
                    b^2/2*(J0b^2 + J1b^2) - b/rho(i)*J0b*J1b -...
                    (a^2/2*(J0a^2 + J1a^2) - a/rho(i)*J0a*J1a));
            else
                J0ai = D(i,1);
                J1ai = E(i,1);
                J0bi = D(i,2);
                J1bi = E(i,2);
                J0aj = D(j,1);
                J1aj = E(j,1);
                J0bj = D(j,2);
                J1bj = E(j,2);
                A(i,j) = C(i)*C(j)*...
                    (-a*rho(i)*J1aj*J0ai...
                    + a*rho(j)*J0aj*J1ai...
                    + b*rho(i)*J1bj*J0bi ...
                    - b*rho(j)*J0bj*J1bi)/(rho(j)^2-rho(i)^2);
            end
            A(j,i) = A(i,j);
            A(j,i) = A(i,j);
        end
        
        % Computing B(1:i-1)*B'(1:i-1)
        if i>1
            BBprime(1:i-2,1:i-2) = BBprime(1:i-2,1:i-2) + B(i-1,1:i-2)*B(i-1,1:i-2)';
            BBprime(:,i-1) = B*B(i-1,:)';
            BBprime(i-1,:) =  B(i-1,:)*B';
        end
        % Next basis vector
        
%         Ei = [ -(B(1:i-1,1:i-1)*B(1:i-1,1:i-1)')*(A(1:i-1,i));1];
        Ei = [ -(BBprime(1:i-1,1:i-1)*A(1:i-1,i));1];
        
        % Norm of this vector
        normSQ(i,1) = sqrt(Ei'*A(1:i,1:i)*Ei);
        Ei = Ei/normSQ(i);
        % Update B
        B(1:i,i) = Ei;
        % Coordinate along Ek
        beta(i,1) = sum(f_ei(1:i).*Ei);
        alpha(1:i,1) = B(1:i,1:i)*beta(1:i);
        resH1square = resH1square - beta(i,1)^2;
        t2 = tic;
        resLinf = norm(storedValsFunc - alpha(1:i,1)'*storedValsJ0Mat(1:i,:),'inf');
        totalLinfNorm = totalLinfNorm + toc(t2);
        if resH1square < tolH1^2
            break
        elseif resLinf < tolInf
            break
        else
            
        end
    end
    totalOrtho = totalOrtho + toc(tOrtho);
    if resH1square < tolH1^2
        
        reachedTol = true;
        P = i;
    elseif resLinf < tolInf
        reachedTol = true;
        P = i;
    else
        P = P+batch_size;
        nextFreqApprox = rho(end) + pi;
    end
    
end
alpha = (B(1:P,1:P)*beta(1:P))./rho(1:P);
w0 = alpha.*C(1:P);
timeTotal = toc(timeTotal);
rho = rho(1:P);
fprintf('Radial quadrature of %d components computed in %d seconds \n',P,timeTotal);
fprintf('Spent %d seconds building the schmidt basis \n',totalOrtho);
fprintf('Spent %d seconds evaluating Bessel functions \n',totalBessel);
fprintf('Spent %d seconds evaluating Linf norm \n',totalLinfNorm);

% Graph
if askGraph
    figure
    quad = @(x)(besselApprox(rho,alpha,0,x,true));
    t = linspace(a/2,(b+1)/2,length(rho)*50);
    tTestLinf = linspace(a,b,length(rho)*50);
    plot(t,log(abs(func(t)-quad(t)+quad(b)-func(b))));
    figure
    plot(rho,alpha);
    
    resLinf = max(abs(quad(tTestLinf) - func(tTestLinf)-quad(b)+func(b)));
    fprintf('Linf error : %d \n',resLinf);
    
end



