loadPackage "Book3264Examples"

Toep = (aa,m,n) -> (
    assert(#aa == m+n-1);
    return matrix for i from 0 to m-1 list for j from 0 to n-1 list aa#(n-1+i-j)
)
Dmatrix = (e,f,k,F) -> ( -- The matrix corresponding to M_k(E->F), where E = O^e, F has rank f
    R := F.AbstractVariety.IntersectionRing;
    cs := chern(f-e+1,f+e-1-2*k,F);
    return Toep(cs,e-k,e-k);
)

recursiveZ = method()
recursiveZ FlagBundle := Gr -> (
    m := Gr.BundleRanks#0 - 1;
    r := Gr.BundleRanks#1;
    R := Gr.IntersectionRing;
    ee := for j from 0 to m+1 list schubertCycle(toList(j:1)|toList((m+1-j):0), Gr);
    h := i -> schubertCycle({i}|toList(m:0), Gr);
    print ee; print h;

    ff := toList((m+1):0_R);
    for i from 0 to r do (
        Z' := sum(1..(m+1),
            j -> (-1)^(j-1) * ee#j * ff#(i-j)
        ) + h(i);
        print("Pass "|toString(i)|": "|toString(Z'));
        ff = insert(i,Z',ff);
    );
    return ff#r;
)
recursiveZ (ZZ,ZZ) := (r,m) -> recursiveZ flagBundle({m+1,r}, VariableNames => {,c})

end

restart
load "Porteous.m2"

(r,m) = (2,2)
Gr = flagBundle({m+1,r}, VariableNames => {,c});
Zrec = recursiveZ(Gr);
--placeholderToSchubertBasis(Zrec, Gr)

e = r+m+2; f = 2*(m+1); k = m+2;
--Gr = flagBundle({m+1,r}, VariableNames => {,c});
(S,Q) = bundles Gr;
F = dual(S+S);
D = Dmatrix(e,f,k,F);
Z = det D;

Z == Zrec

print placeholderToSchubertBasis(Z, Gr);
print placeholderToSchubertBasis(Zrec, Gr);

sigma1 = schubertCycle({1}|toList(m:0), Gr)

placeholderToSchubertBasis(Z * sigma1^r, Gr)
placeholderToSchubertBasis(Zrec * sigma1^r, Gr)