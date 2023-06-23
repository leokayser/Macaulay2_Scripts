load "../Chopnumbers.m2"

n = 2
S = ZZ/2[x,y,z];

findMonIdeals = method()
findMonIdeals (ZZ,ZZ) := (d,numGens) -> (
    i := 0;
    mons := (entries basis(d,S))_0;
    r := binomial(n+d,n) - numGens;
    << ", scheme length = " << r << endl;
    for g in subsets(mons,numGens) do (
        I = monomialIdeal(g);
        if hilbertPolynomial(I, Projective=>false) != r then (
            --<< g << ": wrong HP " << hilbertPolynomial(I, Projective=>false) << endl;
            continue;
        );
        if not verifyESC(saturate(I),r) then (
            << g << ": doesn't satisfy ESC." << endl;
            continue;
        );
        i = i+1;
        << endl;
        << "----> sucess " << g << " <----" << endl;
        << endl;
    );
    print(i);
) 

findMonIdeals(5,3)

--funny := d -> ideal(x^(d+1)*y^d, y^(d+1)*z^d, z^(d+1)*x^d)

-- multMatrixRank(module I, 5, 1)