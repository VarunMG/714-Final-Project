a = -10;
b = 10;
N = 256;
T = 0.5;

tic;
[tVals,xVals,uSpecLie] = KdV_Lie_Godunov(a,b,N,T);
toc;
[tVals,xVals,uSpecStrang] = spectral_Strang_KdV3(a,b,N,T);
toc;
%[tVals,xVals,uSpecSWSS] = spectral_SWSS_KdV(a,b,N,T);
toc;

%% plotting bit

tIters = length(tVals);

ind1 = 1;
ind2 = floor(tIters/3);
ind3 = floor(2*tIters/3);
ind4 = tIters;

figure;

subplot(2,2,1)
hold on;
plot(xVals,uSpecLie(:,ind1))
plot(xVals,uSpecStrang(:,ind1))
%plot(xVals,uSpecSWSS(:,ind1))
title('u at t=0')
xlabel('x')
ylabel('u')
hold off;

subplot(2,2,2)
hold on;
plot(xVals,uSpecLie(:,ind2))
plot(xVals,uSpecStrang(:,ind2))
%plot(xVals,uSpecSWSS(:,ind2))
title('u at t=1/6')
xlabel('x')
ylabel('u')
hold off;

subplot(2,2,3)
hold on;
plot(xVals,uSpecLie(:,ind3))
plot(xVals,uSpecStrang(:,ind3))
%plot(xVals,uSpecSWSS(:,ind3))
title('u at t=1/3')
xlabel('x')
ylabel('u')
hold off;

subplot(2,2,4)
hold on;
plot(xVals,uSpecLie(:,ind4))
plot(xVals,uSpecStrang(:,ind4))
%plot(xVals,uSpecSWSS(:,ind4))
title('u at t=1/2')
xlabel('x')
ylabel('u')
hold off;

legend('Godunov Lie','Spectral Strang')