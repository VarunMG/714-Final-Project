function [tVals,xVals,u] = KdV_Lie(a,b,N,T)
h = 1/N;
dt = 0.01*h;
xVals = linspace(a,b,N);
tVals = 0:dt:T;
[~,tIters] = size(tVals);


r = dt/(2*h);

u = zeros(N,tIters);
u0 = @(x) 1.5+sin(2*pi*x);
u(:,1) = u0(xVals);

for i=1:tIters-1
%     u(2:end,i+1) = u(2:end,i) - 0.5*r*(u(2:end,i).^2-u(1:(end-1),i).^2);
%     u(1,i+1) = u(end,i+1);
    u(2:end-1,i+1) = 0.5*(u(3:end,i)+u(1:(end-2),i)) - 3*r*(u(3:end,i).^2-u(1:(end-2),i).^2);
    u(1,i+1) = 0.5*(u(2,i)+u(end,i)) - 3*r*(u(2,i).^2 - u(end,i).^2);
    u(end,i+1) = 0.5*(u(1,i)+u(end-1,i)) - 3*r*(u(1,i).^2 - u(end-1,i).^2);

%     u(2:end,i+1) = -r*u(2:end,i).*(u(2:end,i)-u(1:(end-1),i))+u(2:end,i);
%     u(1,i+1) = u(end,i+1);
end
end