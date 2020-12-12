A=5;
exact =  @(x,t) 1/2*A*(sech(sqrt(A)/2*(x-A*t+0))).^2;

nVals = [100,150,200,250,300];
errs = zeros(1,5);

for i=1:length(nVals)
    N = nVals(i);
    tic;
    [tVals,xVals,u] = spectral_Lie_KdV2(-10,10,N,0.5);
    toc;
    exactSol = zeros(N,length(tVals));
    for j=1:length(tVals)
        exactSol(:,j) = exact(xVals,tVals(j));
    end  
    err = max(max(abs(exactSol(:,:)-u(:,:))));
    errs(i) = err;
    disp('loop done')
end