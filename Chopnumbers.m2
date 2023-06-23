interestingRange = method(Options => {Igc => true})
interestingRange (ZZ,ZZ) := o -> (n,d) -> (
    rmin := if o.Igc == true then
        ceiling( ((n+1)*binomial(n+d,n) - binomial(n+d+1,n) + 1)/n )
    else
        binomial(n+d-1,n)+1;
    rmax := binomial(n+d,n) - (n+1);
    return (rmin,rmax);
);

regH = (n,r) -> (
    t := 1;
    while binomial(n+t,n) < r do t = t+1;
    return t;
);
idealHF = (n,r,t) -> min(binomial(t+n,n), r);
alternatingSum = (n,r,t) -> (
    d := regH(n,r);
    g := binomial(d+n,n) - idealHF(n,r,d);
    return sum(0..floor(t/d), k -> (-1)^k * binomial(t-k*d+n,n) * binomial(g,k));
);
expectedHF = (n,r,t) -> (
    d = regH(n,r);
    e = expectedGapSize(n,r);
    return if t < d+e then
        alternatingSum(n,r,t)
    else
        idealHF(n,r,t);
)
expectedGapSize = (n,r) -> (
    if r == 1 or (n,r) == (2,4) then return 1;
    d := regH(n,r);
    if r >= binomial(n+d,n)-n then return infinity;
    e := 1;
    while alternatingSum(n,r,d+e) > idealHF(n,r,d+e) do e = e+1;
    return e;
)

igcBeta1 = (n,r) -> (
    d := regH(n,r);
    beta1d := binomial(n+d,n)-r;
    beta1d1 := max(0, binomial(n+d+1,n)-r - beta1d*(n+1));
    return new HashTable from {d=>beta1d, (d+1)=>beta1d1};
);

verifyESC = method() -- This method assumes that I is a 1-dim'l graded ideal of degree r
verifyESC (Ideal,ZZ) := (I,r) -> (
    S := ring(I);
    n := numgens(S)-1;
    d := regH(n,r);
    gapend := d + expectedGapSize(n,r);

    gs := first entries gens I;
    if min(apply(gs, f -> (degree f)_0)) != d then (print("Eq in deg <d."); return false;); -- Check that I has generic HF up to d

    -- The case e=1, only need to check that IGC holds
    if gapend == d+1 then (
        print("IGC case.");
        elapsedTime(mgs := first entries mingens I);
        beta1 = igcBeta1(n,r);
        return number(mgs, f -> (degree f)_0 == d) == beta1#d and number(mgs, f -> (degree f)_0 == d+1) == beta1#(d+1);
    );

    --print("Gappy case.");
    -- The case e>1, need to calculate the Hilbert function by expanding the Hilbert series
    J := ideal select(gs, f -> (degree f)_0 == d);
    --elapsedTime ()
    hs := hilbertSeries(J, Order=>gapend+1);
    T := (ring hs)_0;
    for t from d to gapend do (
        if coefficient(T^t, hs) != expectedHF(n,r,t) then return false;
    );
    return true;
);

-- file = "P3_gaps"
-- for r to 100 do (
--     e = expectedGapSize(3,r);
--     if e == infinity then continue;
--     d = regH(3,r);
--     file << r << " " << d << " " << (d+e) << endl;
-- )