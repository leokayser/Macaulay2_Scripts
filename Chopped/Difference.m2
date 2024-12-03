load "TestSeq.m2"

diffTable = method()
diffTable (ZZ,ZZ) := (n,d) -> (
    rmax := (interestingRange(n,d))_1;

    I := randomPoints(n,rmax);
    --verifyESC (I,rmax);
    g := select(first entries mingens I, f -> (degree f)_0 == d);
    C := ideal(g_(toList(1..n)));

    hI := t -> idealHF(n,rmax,t);
    hC := t -> hilbertFunction({t},C);

    Diff := (f,t) -> f(t)-f(t-1);
    print(hilbertPolynomial(C, Projective=>false));
    l := for t to (d-1)*n+1 list {t, hI(t), hC(t)};
    Dl := for t to (d-1)*n+1 list {t, Diff_hI(t), Diff_hC(t), Diff(k->expectedHF(n,rmax,k),t)};
    return (Dl,l);
);

--Dl = first diffTable(4,5)
n=4; dmax=5;
rmax = last interestingRange(n,dmax);
file = "plotdata/gaps_P" | toString(n);
file << "r d e de" << endl;
for r from 1 to rmax do (
    d = regH(n,r);
    e = expectedGapSize(n,r);
    if e == infinity then continue;
    file << r << " " << d << " " << e << " " << (d+e) << endl;
)
file << close;
