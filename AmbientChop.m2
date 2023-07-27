kk = ZZ/101
needsPackage "NormalToricVarieties"

SampleVariety = new Type of Ideal

makeToricSV = A -> (
    n := numColumns A;
    N := numRows A;
    Pn := toricProjectiveSpace(n, CoefficientRing => kk);
    PN := toricProjectiveSpace(N, CoefficientRing => kk);
    IA := ideal map(PN, Pn, A);
    monMap := mat -> matrix{{1} | for i to numRows(A)-1 list product(flatten entries mat,flatten entries A^{i}, (x,a)->x^a)};
    params = new HashTable from {param => monMap, dim => n, codim => N-n};
    sv := new SampleVariety from merge(IA, params, );
    return sv;
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
productMatrix(pn(2),veroneseMatrix(1,2))

sv = makeToricSV veroneseMatrix(2,2)

rat = d -> matrix(for i from 1 to d list {i})
makeToricSV rat 1

secantLin = method()
secantLin (Matrix,Ring) := (mat,S) -> (
    minorSize := min(numRows(mat)+1, numColumns(mat));
    return ideal mingens minors(minorSize, sub(mat,S) || vars(S));
)
(n,d) = (2,4);

vd = veronese(n,d,K); IX = ker(vd); S = ring IX;
rmax = codim(IX) -- So that Span(Z) \cap X = Z (generically) by Multisecant lemma;

r = 9;
Z = mapPoints(vd,random(K^r,K^(n+1)))

L = secantLin(Z,S);
J = IX + L;

tmax = 2;
netList transpose reverse for t to tmax list {tmax-t, hilbertFunction(tmax-t,J)} -- Calling hilbertFunction in reverse to avoid redundant computation
