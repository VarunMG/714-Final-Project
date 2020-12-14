A=5;
exact =  @(x,t) 1/2*A*(sech(sqrt(A)/2*(x-A*t+0))).^2;

%nVals = [120,160,200,240,280,320];
nVals = [100,150,200,250,300];
%nVals = [100,125,150,200,250];

iters = length(nVals);

errsSpecLie = zeros(1,iters);
errsSpecStrang = zeros(1,iters);
errsSpecSWSS = zeros(1,iters);

massDiffSpecLie = zeros(1,iters);
massDiffSpecStrang = zeros(1,iters);
massDiffSpecSWSS = zeros(1,iters);

momDiffSpecLie = zeros(1,iters);
momDiffSpecStrang = zeros(1,iters);
momDiffSpecSWSS = zeros(1,iters);

for i=1:iters
    N = nVals(i);
    tic;
    [tVals,xVals,uSpecLie] = spectral_Lie_KdV2(-10,10,N,0.5);
    [tVals,xVals,uSpecStrang] = spectral_Strang_KdV3(-10,10,N,0.5);
    [tVals,xVals,uSpecSWSS] = spectral_SWSS_KdV(-10,10,N,0.5);
    toc;
    [X,T] = meshgrid(xVals,tVals);
    exactSol = exact(X,T); 
    
    errSpecLie = max(max(abs(exactSol(:,:)-uSpecLie(:,:).')));
    errSpecStrang = max(max(abs(exactSol(:,:)-uSpecStrang(:,:).')));
    errSpecSWSS = max(max(abs(exactSol(:,:)-uSpecSWSS(:,:).')));
    
    mSpecLie = mass(uSpecLie);
    mSpecStrang = mass(uSpecStrang);
    mSpecSWSS = mass(uSpecSWSS);
    
    pSpecLie = momentum(uSpecLie);
    pSpecStrang = momentum(uSpecStrang);
    pSpecSWSS = momentum(uSpecSWSS);
    
    mDiffSpecLie = abs(max(mSpecLie)-min(mSpecLie));
    mDiffSpecStrang = abs(max(mSpecStrang)-min(mSpecStrang));
    mDiffSpecSWSS = abs(max(mSpecSWSS)-min(mSpecSWSS));
    
    pDiffSpecLie = abs(max(pSpecLie)-min(pSpecLie));
    pDiffSpecStrang = abs(max(pSpecStrang)-min(pSpecStrang));
    pDiffSpecSWSS = abs(max(pSpecSWSS)-min(pSpecSWSS));
    

    massDiffSpecLie(i) = mDiffSpecLie;
    massDiffSpecStrang(i) = mDiffSpecStrang;
    massDiffSpecSWSS(i) = mDiffSpecSWSS;
    
    momDiffSpecLie(i) = pDiffSpecLie;
    momDiffSpecStrang(i) = pDiffSpecStrang;
    momDiffSpecSWSS(i) = pDiffSpecSWSS;
    
    errsSpecLie(i) = errSpecLie;
    errsSpecStrang(i) = errSpecStrang;
    errsSpecSWSS(i) = errSpecSWSS;
    disp('loop done')
end
%% plotter stuff
figure;
hold on;
plot(log(nVals),log(errsSpecLie));
plot(log(nVals),log(errsSpecStrang));
plot(log(nVals),log(errsSpecSWSS));
hold off;
legend('Spectral Lie','Spectral Strang','Spectral SWSS')
title(' log(error) vs. log(N)')
xlabel('log(N)');
ylabel('log(error)')

disp('log(N) vs log(error) linear fit')
polyfit(log(nVals),log(errs),1)

figure;
plot(log(nVals),log(massDiffs))
title('log of Mass Variation vs. log(N)')
xlabel('log(N)')
ylabel('log(Mass Variation)')
disp('log(Mass Variation) vs. log(N) linear fit')
polyfit(log(nVals),log(massDiffs),1)

figure;
plot(log(nVals),log(momentumDiffs))
title('log(Momentum Variation) vs. log(N)')
xlabel('log(N)')
ylabel('log(Momentum Variation)')
disp('log(Momentum Variation) vs. log(N) linear fit')
polyfit(log(nVals),log(momentumDiffs),1)