restart
loadPackage "SymmetricPolynomials"
loadPackage "SpechtModule"

m = 2; R = QQ[x_0..x_m]
R' = QQ[x_0..x_m,y_0..y_m]
phi = map(R,R',gens R | gens R)
h = (j,R) -> sum flatten entries basis(j,R)
elementarySymmetric phi h(2,R')

ee = {1_R} | elementarySymmetricPolynomials R 

hhdecomp = r -> (
    ff := toList((m+1):0_R);
    for i from 0 to r do (
        Z' := sum(1..(m+1), j -> (-1)^(j-1) * ee_j * ff#(i-j)) + h(i,R);
        ff = insert(i,Z',ff);
    );
    return ff#r;
)
hhdecomp 2 == phi h(2,R')