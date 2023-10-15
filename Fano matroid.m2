loadPackage "Matroids"
loadPackage "NormalToricVarieties"

A = sub(transpose matrix{{1,0,0},{1,0,1},{0,1,0},{0,0,1},{1,1,1},{0,1,1},{1,1,0}},ZZ/2)
M = matroid A 
C = circuits M
B = bases M
M'= dual M
C'= circuits M'

S = QQ[x_0..x_6]
f = sum(B, b -> product(toList b, i -> x_i))
f' = flatten entries diff(vars S, f)
Hs = apply(f', g -> diff(vars S,diff(transpose vars S, g)))
lambdas = eigenvalues sub(Hs_0,QQ)

V = basisIndicatorMatrix M
X = normalToricVariety V