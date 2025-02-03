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
end

restart
load "Porteous.m2"

m = 1
results = for r from 1 to 10 list (
    -- O^{r+m+2} -> (S^*)^2
    e = r+m+2; f = 2*(m+1); k = m+2; (e,f,k);
    Gr = flagBundle({m+1,r}, VariableNames => {,c});
    (S,Q) = bundles Gr;
    F = dual(S+S);
    D = Dmatrix(e,f,k,F);
    {r, placeholderToSchubertBasis(det D, Gr)}
);
netList results


