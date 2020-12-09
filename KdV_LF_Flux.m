function [tVals,xVals,u] = KdV_LF_Flux(a,b,T,n)
%n is number of spatial grid points
%m is number of time steps
%solving on x \in [a,b] and t in [0,T]
u0 = @(x) 1.5+ sin(2*pi*x);
xVals = linspace(a,b,n);
dx = xVals(2)-xVals(1);
dt = 0.01*dx;
tVals = 0:dt:T;
[~,tIters] = size(tVals);
r = dt/(2*dx);
u = zeros(n,tIters);
u(:,1) = u0(xVals);
for i=1:(tIters-1)
    u(3:(end-1),i+1) = u(3:(end-1),i) + r*( 3*u((3:end-1),i).^2 - ...
        3*(u(1:(end-3),i)).^2) - 0.5*(u(3:(end-1),i) - u(2:(end-2),i));
        %- (1/(dx*dx*dx))*(u(4:end,i) - 2*u(3:(end-1),i) + 2*u(2:(end-2),i) - u(1:(end-3),i));
    u(1,i+1) = u(1,i) + r*(3*u(1,i)^2-3*u(end-1,i)^2) - 0.5*(u(1,i) - u(end-1,i));
        %-(1/(dx*dx*dx))*(u(2,i)-2*u(1,i)+2*u(end,i)-u(end-1,i));
    u(2,i+1) = u(2,i) + r*(3*u(2,i)^2-3*u(end,i)^2) - 0.5*(u(2,i) - u(end,i));
        %-(1/(dx*dx*dx))*(u(3,i)-2*u(2,i)+2*u(1,i)-u(end,i));
    u(end,i+1) = u(end,i) + r*(3*u(end,i)^2-3*u(end-2,i)^2) - 0.5*(u(end,i)-u(end-2,i));
        %-(1/(dx*dx*dx))*(u(1,i)-2*u(end,i)+2*u(end-1,i)-u(end-2,i));
end
%maxU = max(u(:));
%minU = min(u(:));
for j=1:(tIters-1)
    plot(xVals,u(:,j))
    ylim([0,3]);
    pause(0.05)
end
end