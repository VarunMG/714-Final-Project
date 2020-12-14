function [tVals,xVals,u] = KdV_Lie_LF(a,b,N,T)
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

u0 = @(x) f(x);

%u0 = @(x) sinc(2*x);

u = zeros(N,tIters);
u(:,1) = u0(xVals);
uHat = zeros(N,tIters);
uHat(:,1) = fft(u(:,1));


ik3 = 1j*k.^3;
eik3dt = exp(ik3*dt);

r = dt/(2*dx);

for n=1:(tIters-1)
    u1 = zeros(N,1);
    u1(2:end-1) = 0.5*(u(3:end,n)+u(1:(end-2),n)) - 3*r*(u(3:end,n).^2-u(1:(end-2),n).^2);
    u1(1) = 0.5*(u(2,n)+u(end,n)) - 3*r*(u(2,n).^2 - u(end,n).^2);
    u1(end) = 0.5*(u(1,n)+u(end-1,n)) - 3*r*(u(1,n).^2 - u(end-1,n).^2);
    
    u(:,n+1) =  real(ifft(fft(u1).*eik3dt.'));
end
end


% function [tVals,xVals,u] = KdV_Lie_LF(a,b,N,T)
% dx = 1/N;
% dt = 0.1/(N^2);
% xVals = linspace(a,b,N);
% tVals = 0:dt:T;
% [~,tIters] = size(tVals);
% dk = 2*pi/(N*dx);
% k = [0:dk:N/2*dk,-(N/2-1)*dk:dk:-dk];
% 
% r = dt/(2*dx);
% 
% u = zeros(N,tIters);
% 
% A=5;
% f =  @(x) 1/2*A*(sech(sqrt(A)/2*(x+0))).^2;
% 
% B=20;
% g = @(x) 1/2*B*(sech(sqrt(B)/2*(x+4))).^2;
% 
% u0 = @(x) f(x);
% 
% u(:,1) = u0(xVals);
% 
% ik3 = 1j*k.^3;
% eik3dt = exp(ik3*dt);
% 
% for n=1:tIters-1
% %     u(2:end,i+1) = u(2:end,i) - 0.5*r*(u(2:end,i).^2-u(1:(end-1),i).^2);
% %     u(1,i+1) = u(end,i+1);
%     u1 = zeros(N,1);
%     u1(2:end-1) = 0.5*(u(3:end,n)+u(1:(end-2),n)) - 3*r*(u(3:end,n).^2-u(1:(end-2),n).^2);
%     u1(1) = 0.5*(u(2,n)+u(end,n)) - 3*r*(u(2,n).^2 - u(end,n).^2);
%     u1(end) = 0.5*(u(1,n)+u(end-1,n)) - 3*r*(u(1,n).^2 - u(end-1,n).^2);
%     
%     u(:,n+1) = real(ifft(fft(u1).*eik3dt.'));
% 
% %     u(2:end,i+1) = -r*u(2:end,i).*(u(2:end,i)-u(1:(end-1),i))+u(2:end,i);
% %     u(1,i+1) = u(end,i+1);
% end
% end