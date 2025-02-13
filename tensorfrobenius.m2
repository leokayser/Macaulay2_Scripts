restart
p = 2;
(d1,d2) = (2,3);
d = d1+d2-2
KK = frac (ZZ/p[a_0..a_d])
S = KK[x,y]
F = sum(0..d, i -> a_i*x^(d-i)*y^i)
I = inverseSystem(F)