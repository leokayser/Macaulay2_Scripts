interestingRange = method()
interestingRange (ZZ,ZZ) := (n,d) -> (
    rmin := ceiling( ((n+1)*binomial(n+d,n) - binomial(n+d+1,n) + 1)/n );
    rmax := binomial(n+d,n) - (n+1);
    return (rmin,rmax);
)

idealHF = (n,r,j) -> binomial(j+n,n) - r;
alternatingSum = (n,d,r,e) -> sum(floor(e/d)+1, k -> (-1)^k * binomial(e-k*d+n,n) * binomial(idealHF(n,r,d),k+1));
expectedHF = (n,d,r,e) -> min(alternatingSum(n,d,r,e), idealHF(n,r,d+e));

expectedGapSize = (n,d,r) -> (
    e := 1;
    while alternatingSum(n,d,r,e) < idealHF(n,r,d+e) do e = e+1;
    return e;
)

chopDegree = (n,r) -> (
    j := 1;
    while binomial(n+j,n) <= r do j = j+1;
    return j;
)

testCase = method()
testCase (PolynomialRing,ZZ) := (S,r) -> (
    kk := coefficientRing S;
    n := (numgens S)-1;
    d := chopDegree(n,r);

    identStr := concatenate("(d=", toString(d), ",r=", toString(r), ")");

    e := expectedGapSize(n,d,r);
    print(identStr | " -> expected gap = " | e);
    
    while true do (
        randPts := for i from 1 to r list matrix{for j to n list random(kk)};
        I := intersect(apply(randPts, pt -> minors(2, pt || basis(1,S))));
        if #((entries super basis(d,I))_0) != idealHF(n,r,d) then continue;

        actualHF := e -> (
            SeId := ideal(basis(e,S))*I;
            return #((entries super basis(d+e, SeId))_0);
        );
        if actualHF(e) == idealHF(n,r,d+e) and  actualHF(e-1) == alternatingSum(n,d,r,e-1) then (
            print(identStr | " valid!");
            return randPts;
        );
        print(identStr | " invalid, trying again");
    );
)
