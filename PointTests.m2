restart
notify = true
load "PointIdeals.m2"

r = 25
d = (expectedNumGens(r))_0_0
Z = randomQPoints(r)
I = IGenD(Z,d);


M = R^1/I;
Mred = R^1/vanishingIdeal(Z);

for i to 12 list hilbertFunction(i, I)
F = res M
F' = res Mred
betti F
betti F'
A = F'.dd_2;
netList transpose apply(numcols A, i -> (apply(numrows A, j -> degree A_i_j)))

netList for r from 1 to 50 list (
    g := expectedNumGens(r);
    if not isInteresting(r) then continue;
    e := g_0_0 + 1;
    while expectedHilbFun(r,g,e) > r do e = e+1;
    {r} | g | {(0..e) / (j->expectedHilbFun(r,j)) | (1:"...")}
)


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