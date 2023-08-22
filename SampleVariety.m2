needsPackage "NormalToricVarieties"
needsPackage "Points"

SampleVariety = new Type of Ideal
-- Keys of a SampleVariety (in addition to those of an Ideal):
-- sample:   a function which returns a random kk-point on the variety
-- codim:    the codimension of X (we assume that X is linearly non-degenerate)

sampleRandomPoint = method(Options => {MaxTrials => 3})
sampleRandomPoint (SampleVariety) := o -> sv -> (
    kk := coefficientRing ring sv;
    for t from 1 to o.MaxTrials do (
        try (
            p := random(kk^1,kk^(sv.paramdim));
            z := sv.param(p)
        ) then (
            return z;
        ) else (
            continue;
        );
    );
    throw error "Failed to produce a kk-point.";
)

-- AFFINE Graßmannian convention (!) grassSV(k,n) -> Gr(k,n) = PGr(k-1,n-1)
grassSV = method(Options => {kk => ZZ/101})
grassSV (ZZ,ZZ) := o -> (k,n) -> (
    I := Grassmannian(k-1, n-1, CoefficientRing=>o.kk, Variable=>x);
    pluckerEmb := () -> (
        A := id_(o.kk^k) | random(o.kk^k, o.kk^(n-k));
        return matrix{apply(subsets(n,k), s -> det submatrix(A, s))};
    );
    params := new HashTable from {sample => pluckerEmb, codim => binomial(n,k)-1-k*(n-k)};
    return new SampleVariety from merge(I, params, );
)

toricSV = method(Options => {kk => ZZ/101}) -- Requires an integer matrix in "affine convention"
toricSV (Matrix) := o -> A -> (
    n := numColumns A;
    N := numRows A;
    An := affineSpace(n, CoefficientRing => o.kk);
    PN := toricProjectiveSpace(N, CoefficientRing => o.kk);
    IA := ideal map(PN, An, A);
    monMap := () -> matrix{{1} | for i to numRows(A)-1 list product(flatten entries random(o.kk^1,o.kk^n), flatten entries A^{i}, (x,a)->x^a)};
    params := new HashTable from {sample => monMap, codim => N-n};
    return new SampleVariety from merge(IA, params, );
)

rankSV = method(Options => {kk => ZZ/101})
rankSV (ZZ,ZZ,ZZ) := o -> (n1,n2,r) -> (
    S := o.kk[x_(0,0)..x_(n1-1,n2-1)];
    Ir := minors(r+1, genericMatrix(S,n1,n2));
    tensorSum := () -> flatten sum(1..r, i -> random(o.kk^n1,o.kk^1) * random(o.kk^1,o.kk^n2));
    params := new HashTable from {sample => tensorSum, codim => n1*n2-r*(n1+n2-r)};
    return new SampleVariety from merge(Ir, params, );
)

affinizeMatrix = A -> submatrix'(A,{0},{0})
projectivizeMatrix = A -> (
    a := entries A;
    d := max(a/sum);
    firstRow := {d} | toList(#a_0:0);
    return matrix( {firstRow} | apply(a, row -> {d - sum(row)} | row));
)

veroneseMatrix = (n,d) -> affinizeMatrix matrix( flatten entries super basis(d,ZZ[vars(0..n)]) / exponents / flatten )
pn = n -> veroneseMatrix(n,1)
productMatrix = M -> (
    if instance(M,Matrix) then return M;
    if #M == 1 then return M_0;
    m := toList(M) / projectivizeMatrix / entries;
    p := fold((a,b) -> (a ** b) / toList / flatten, m);
    colsToRemove = {0} | (accumulate(plus, 0, m/first/length));
    return submatrix'(matrix p, {0}, colsToRemove);
)

nilpotentSV = method(Options => {kk => ZZ/101})
nilpotentSV (ZZ) := o -> n -> (
    S := o.kk[x_(0,0)..x_(n-1,n-1)];
    use S[t];
    charpol := determinant(t*id_(S^n) - genericMatrix(S,n,n));
    coeffs := sub(last coefficients(charpol, Monomials=>for i to n-1 list t^i), S);
    I := ideal coeffs;
    -- We need to project to a polynomial ring in n^2-1 variables to make V(I) non-degenerate.
    S' := o.kk[drop(gens S,-1)];
    projection := map(S', S, gens S'|{-sum(n-1, i->x_(i,i))});
    I' := projection I;

    Jn := matrix for i from 1 to n list for j from 1 to n list if j==i+1 then 1 else 0;
    conjMap := mat -> (
        A := reshape(o.kk^n, o.kk^n, mat);
        N := A * Jn * A^-1;
        return submatrix'(flatten N, {n^2-1});
    );
    params := new HashTable from {param => conjMap, paramdim => n^2, codim => n-1};
    return new SampleVariety from merge(I', params, );
)

linSpan = method()
linSpan (List,Ring) := (Z,S) -> (
    mat := fold(Z, (z1,z2) -> z1||z2);
    minorSize := min(numRows(mat)+1, numColumns(mat));
    mat' := sub(mat,S) || vars(S);
    return ideal mingens minors(minorSize, mat');
)

toPoints = (Z,S) -> transpose sub(fold(Z, (z,z') ->  z||z'), S)
genPtsIdeal = (S,r) -> intersect for i from 0 to r-1 list ideal drop(gens S, {i,i})


hilbFunTable = args -> (
    tmax := last args;
    Ms := toList drop(args, -1);
    return transpose for t to tmax list {t} | apply(Ms, M -> HF(M,t));
)

disp = (sv,tmax) -> (
    r := sv#codim;
    Z := for i from 1 to r list sampleRandomPoint(sv);

    L := linSpan(Z,ring sv);
    J := sv + L;
    IZ := points toPoints(Z, ring sv);

    print netList apply({{""},{" S "},{" I "},{" L "},{"I+L"},{"I_Z"}}, hilbFunTable(module ring sv, module sv, module L,module J, module IZ, tmax), (a,b)->a|b);
    return J;
)

HF = method()
HF (Module,ZZ) := (M,t) -> numColumns super basis(t,M)

