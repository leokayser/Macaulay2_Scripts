restart


multMatrixRank = method()
multMatrixRank (Module,ZZ,ZZ) := (M, d, e) -> (
    S := ring(M);
    SM := ideal(basis(e,S))*M;
    return #(entries super basis(d+e, SM))_0
)

n = 2
S = ZZ/2[x,y,z];

findMonIdeals = method()
findMonIdeals (ZZ,ZZ,ZZ) := (d,numGens,e) -> (
    i := 0;
    mons := (entries basis(d,S))_0;
    schemeLength := binomial(n+d,n) - numGens;
    maxRank := min(numGens*binomial(n+e,n), binomial(n+d+e,n) - schemeLength);
    << "maxRank = " << maxRank << ", scheme length = " << schemeLength << endl;
    for g in subsets(mons,numGens) do (
        I = ideal(g);
        if multMatrixRank(module I, d, e) < maxRank or hilbertPolynomial(I, Projective=>false) != schemeLength then
            continue;
        i = i+1;
        print(I);
        print(hilbertPolynomial(I,Projective => false));
    );
    print(i);
) 

findMonIdeals(5,3,2)

funny := d -> ideal(x^(d+1)*y^d, y^(d+1)*z^d, z^(d+1)*x^d)

-- multMatrixRank(module I, 5, 1)