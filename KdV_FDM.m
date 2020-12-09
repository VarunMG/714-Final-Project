function [xVals,tVals,u] = KdV_FDM(a,b,T,n)
xVals = linspace(a,b,n);
dx = xVals(2)-xVals(1);
xVals = [a-2*dx, a-2*dx, xVals, b+dx, b+2*dx];
dt = 0.1*dx*dx*dx;
tVals = 0:dt:T;
[~,tIters] = size(tVals);
u = zeros(n+4,tIters);
u0 = @(x) x.*(1-x);
u(3:(end-2),1) = u0(xVals(3:(end-2)));
for m=1:(tIters-1)
    nonlin = (1/dx).*u(3:(end-2),m).*(u(4:(end-1),m)-u(3:(end-2),m));
    lin = - (1/(dx^3)).*(u(5:end,m)-2*u(4:(end-1),m)+2*u(2:(end-3),m)-u(1:(end-4),m));
    u(3:(end-2),m+1) = nonlin + lin;
end
for j=1:(tIters-1)
    plot(xVals,u(:,j))
    pause(0.05)
end
end