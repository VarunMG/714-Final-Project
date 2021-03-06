function [tVals,xVals,u] = spectral_Strang_KdV3(a,b,N,T)
%%spectral code using a Strang spitting in solving the PDE after transforming
%%into Fourier space 
dt = 0.4/N^2;
tVals = 0:dt:T;
xVals = linspace(a,b,N);
dx = xVals(2)-xVals(1);
[~,tIters] = size(tVals);
dk = 2*pi/(N*dx);
k = [0:dk:N/2*dk,-(N/2-1)*dk:dk:-dk];

%two soliton initial data
A=30;
f =  @(x) 1/2*A*(sech(sqrt(A)/2*(x+8))).^2;

B=20;
g = @(x) 1/2*B*(sech(sqrt(B)/2*(x+4))).^2;

u0 = @(x) f(x) + g(x);


u = zeros(N,tIters);
u(:,1) = u0(xVals);
uHat = zeros(N,tIters);
uHat(:,1) = fft(u(:,1));


ik3 = 1j*k.^3;
eik3dt = exp(ik3*dt);

for n=1:(tIters-1)
    u1 = uHat(:,n) - 0.5*dt*(3i*k.'.*fft(real(ifft(uHat(:,n))).^2));
    u2 =  u1.*eik3dt.';
    uHat(:,n+1) =  u2 - 0.5*dt*(3i*k.'.*fft(real(ifft(u2)).^2));
    u(:,n+1) = real(ifft(uHat(:,n+1)));
end
end