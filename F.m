function val = F(Qi,Qi1)
f = @(x) 3*x.^2;
s = (f(Qi)-f(Qi1))/(Qi-Qi1);
qs = 0;
if Qi1 > qs && s > 0
    val = f(Qi1);
elseif Qi < qs && s <0
    val = f(Qi);
else
    val = f(qs);
end