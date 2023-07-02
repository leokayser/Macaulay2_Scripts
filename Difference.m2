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

Dl = first diffTable(4,5)

-- file = "DiffPlot45.dat"
-- file << "t zero hI hC" << endl;
-- for tup in Dl do (
--     file << tup_0 << " 0 " << tup_1 << " " << tup_2 << endl;
-- )
-- file << close;
