function m = mass(u)
[~,tIters] = size(u);
m = zeros(tIters,1);
for i=1:tIters
    m(i) = sum(u(:,i));
end
end