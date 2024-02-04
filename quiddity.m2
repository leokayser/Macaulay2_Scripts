eta = c -> matrix{{c,-1},{1,0}}

quIdeal = (lambda,m) -> (
    baseRing := if instance(lambda,ZZ) then QQ else ring(lambda);
    R := baseRing[vars(0..(m-1))];
    A := product(gens R, c->eta(c)) - matrix{{lambda,0},{0,lambda}};
    return ideal flatten entries A;
)

glProdIdeal = m -> (
    R := QQ[vars(0..(4*m-1))];
    A := product(0..(m-1), k -> genericMatrix(R,(gens R)_(4*k),2,2)) - matrix{{1,0},{0,1}};
    return ideal flatten entries A;
)

slProdIdeal = m -> (
    I := glProdIdeal m;
    R := ring I;
    sl := ideal apply(0..(m-1), k -> det(genericMatrix(R,(gens R)_(4*k),2,2))-1);
    return I + sl;
)

quBIdeal = m -> (
    R := QQ[vars(0..(m-1))];
    P := product(gens R, c->eta(c));
    return ideal(P_(1,0),P_(0,1),P_(0,0)-P_(1,1));
)

quGrIdeal = (lambda,m) -> (
    baseRing := if instance(lambda,ZZ) then QQ else ring(lambda);
    R := baseRing[vars(-1..(m-1))];
    A := product((gens R)_{1..m}, c->eta(c)) - matrix{{lambda,0},{0,lambda}};
    return homogenize(ideal flatten entries A, (gens R)_0);
)

quBGrIdeal = m -> (
    R := QQ[vars(-1..(m-1))];
    P := product((gens R)_{1..m}, c->eta(c));
    return ideal(homogenize(P_(1,0),P_(0,1),P_(0,0)-P_(1,1)), (gens R)_0);
)

--loadPackage "Depth"

for m from 3 to 20 do (
    --if m%2==1 then continue;
    I := quBGrIdeal m;
    R := ring I;
    print concatenate("m = ", toString m , ", deg = ", toString degree I, ", codim = ", toString codim I);
    --if not isCM(R/I) then print("Not CM?!");
    print minimalBetti I;
    print "";
)


loadPackage "Divisor"

-- https://oeis.org/A006918
for m from 3 to 20 do (
    I := quIdeal(-1,m);
    print(m, degree I);
)
-- -1: {1, 1, 2, 5, 8, 14, 20, 30, 40, 55, 70, 91, 112, 140}
--  1: {0, 1, 2, 5, 8, 14, 20, 30, 40, 55, 70, 91, 112, 140}

-- -1: They are all normal (up to say m=15).
--     m=6 seems to be the only singular case so far...


-- https://oeis.org/A006584
use QQ[lambda]
for m from 2 to 15 list isPrime(quIdeal(-1,m))
--     {1, 2, 4, 10, 16, 28, 40, 60, 80, 110, 140, 182, 224, 280} this is the sum of the above two sequences
