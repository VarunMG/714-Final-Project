function [tVals,xVals,u] = spectral_KdV(N,T)
%Make N a power of 2 or worst case EVEN
%solves on [-pi,pi]
dx = 1/(N-1);
dt = 0.01/(N^2);
tVals = 0:dt:T;
xVals = linspace(-10,10,N);
[~,tIters] = size(tVals);
dk = 2*pi/(N*dx);
k = [0:dk:N/2*dk,-(N/2-1)*dk:dk:-dk];

%initial data
A=20;
u0 =  @(x) 1/2*A*(sech(sqrt(A)/2*(x+8))).^2;


u = zeros(N,tIters);
u(:,1) = u0(xVals);
uHat = zeros(N,tIters);
uHat(:,1) = fft(u(:,1));

ik3 = 1i*k.^3;

for n=1:(tIters-1)
    t = dt*n;
    emik3t = exp(-ik3*dt);
    eik3t = exp(ik3*dt);
    uHat(:,n+1) = uHat(:,n+1) - 0.5*dt*1i*emik3t.'.*k.'.*fft((real(ifft(eik3t*uHat(:,n))).^2));
    u(:,n+1) = real(ifft(uHat(:,n+1)));
end
end