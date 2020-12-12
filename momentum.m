function p = momentum(u)
[~,tIters] = size(u);
p = zeros(tIters,1);
for i=1:tIters
    p(i) = sum(u(:,i).^2);
end
end