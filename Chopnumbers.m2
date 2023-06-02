interestingRange = method()
interestingRange (ZZ,ZZ) := (n,d) -> (
    rmin := ceiling( ((n+1)*binomial(n+d,n) - binomial(n+d+1,n) + 1)/n );
    rmax := binomial(n+d,n) - (n+1);
    return (rmin,rmax);
)

idealHF = (n,r,j) -> min(binomial(j+n,n), r);
alternatingSum = (n,d,r,j) -> (
    g := binomial(d+n,n) - idealHF(n,r,d);
    return sum(0..floor(j/d), k -> (-1)^k * binomial(j-k*d+n,n) * binomial(g,k));
);
expectedHF = (n,d,r,j) -> max(alternatingSum(n,d,r,j), idealHF(n,r,j));

regH = (n,r) -> (
    j := 1;
    while binomial(n+j,n) < r do j = j+1;
    return j;
)

expectedGapSize = (n,r) -> (
    d := regH(n,r);
    if r >= binomial(n+d,n)-n then return infinity;
    e := 1;
    while alternatingSum(n,d,r,d+e) > idealHF(n,r,d+e) do e = e+1;
    return e;
)

randomPoints = method()
randomPoints (PolynomialRing,ZZ) := (S,r) -> (
    n := numgens(S)-1;
    kk := coefficientRing(S);
    return for i from 1 to r list matrix{for j to n list random(kk)};
)

vanishIdeal = method()
vanishIdeal (PolynomialRing,List) := (S,pts) -> intersect(apply(pts, pt -> minors(2, pt || basis(1,S))))

verifyConj = method()
verifyConj (Ideal,ZZ) := (I,r) -> (
    S := ring(I);
    n := numgens(S)-1;
    d := regH(n,r);
    J := ideal super basis(d,I);
    gapend := d + expectedGapSize(n,r);
    hs := hilbertSeries(J, Order=>gapend+1);
    T := (ring hs)_0;
    for j from d to gapend do (
        if coefficient(T^j, hs) > expectedHF(n,d,r,j) then return false;
    );
    return true;
)
