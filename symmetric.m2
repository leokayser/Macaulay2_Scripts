restart
loadPackage "SymmetricPolynomials"
loadPackage "SpechtModule"

m = 2; R = QQ[x_0..x_m]
R' = QQ[x_0..x_m,y_0..y_m]
phi = map(R,R',gens R | gens R)
h = (j,R) -> sum flatten entries basis(j,R)
elementarySymmetric phi h(2,R')

ee = {1_R} | elementarySymmetricPolynomials R 

-- hhdecomp = r -> (
--     ff := toList((m+1):0_R);
--     for i from 0 to r do (
--         Z' := sum(1..(m+1), j -> (-1)^(j-1) * ee_j * ff#(i-j)) + h(i,R);
--         ff = insert(i,Z',ff);
--     );
--     return ff#r;
-- )
-- hhdecomp 3 == phi h(3,R')

h' = i -> (
    if i < 0 then return 0;
    if i == 0 then return 1;
    return sum(1..(m+1), j -> (-1)^(j-1) * ee#j * h'(i-j));
)

for i from 0 to 5 do (print h'(i))