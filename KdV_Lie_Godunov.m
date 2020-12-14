function [tVals,xVals,u] = KdV_Lie_Godunov(a,b,N,T)
%%spectral code using a Lie spitting in solving the PDE after transforming
%%into Fourier space 
dt = 0.1/N^2;
tVals = 0:dt:T;
xVals = linspace(a,b,N);
dx = xVals(2)-xVals(1);
[~,tIters] = size(tVals);
dk = 2*pi/(N*dx);
k = [0:dk:N/2*dk,-(N/2-1)*dk:dk:-dk];


%two soliton initial data
A=5;
f =  @(x) 1/2*A*(sech(sqrt(A)/2*(x+0))).^2;

B=20;
g = @(x) 1/2*B*(sech(sqrt(B)/2*(x+4))).^2;

u0 = @(x) f(x) + g(x);

%u0 = @(x) sinc(2*x);

u = zeros(N,tIters);
u(:,1) = u0(xVals);
uHat = zeros(N,tIters);
uHat(:,1) = fft(u(:,1));


ik3 = 1j*k.^3;
eik3dt = exp(ik3*dt);

r = dt/(dx);

for n=1:(tIters-1)
    u1 = zeros(N,1);
    u1(1) = u(1,n) - r*(F(u(2,n),u(1,n))-F(u(1,n),u(end,n)));
    for j=2:N-1
        u1(j)= u(j,n) - r*(F(u(j+1,n),u(j,n))-F(u(j,n),u(j-1,n)));
    end
    u1(end) = u(end,n) - r*(F(u(1,n),u(end,n))-F(u(end,n),u(end-1,n)));
    
    u(:,n+1) =  real(ifft(fft(u1).*eik3dt.'));
end
end