restart
notify = true
load "HilbertTest.m2"

Z = randomQPoints(3)
vanishingIdeal(Z)

d = 7
s_d = binomial(n+d,n)
r = 33
Z = randomQPoints(r)
I = IgenD(Z,d);

M = R^1/I;
Mred = R^1/Jall(Z);
for i to 14 list hilbertFunction(i, M)
F = res M

betti F


dr = method()
dr := r -> (d := 1; while binomial(n+d,n) < r do (d = d+1); return d)


isInteresting(9)

netList for i from 1 to 100 list (ds = expectedNumGens(i); if not isInteresting(ds) then continue; {i,ds})


S = QQ[z,x,y]
e = {1,2,2}
f = {2,2,2}
ms = for i from 0 to 3 list product(toList apply(0..(i-1),i -> x^(e_i)))*product(toList apply(i..2,i -> y^(f_i)))
A = S^1/(ideal(ms))
F = res A
I = ideal(ms)
netList for i from 5 to 7 list {i, flatten entries super basis(i,I)}
flatten entries super basis(6,I)

(exponents ms_1)_0_1

fm = method()
fm := m -> eta -> product(1..2, i -> product(0..exponents(m)_0_i, x_i - ))