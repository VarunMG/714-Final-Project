function val = visualize_KdV(tVals,xVals,u)
[~,tIters] = size(tVals);
yMax = max(max(abs(u(:,:))));
figure;
for i=1:100:tIters
    plot(xVals,u(:,i))
    ylim([-yMax,yMax]);
    pause(0.01);
end
end