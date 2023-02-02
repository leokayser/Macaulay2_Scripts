n = 2
R = QQ[x_0..x_n]

randomQPoints = method(Options => true)
randomQPoints ZZ := {UseStdPoints => true, IntRange => (-42,42)} >> o -> r -> (
    if o.UseStdPoints then (
        stdPts = apply(entries gens QQ^(n+1), e -> matrix{e});
        stdPts = append(stdPts, sum(stdPts));
        if r <= n+2 then return take(stdPts,r);
    ) else stdPts = {};
    randPts = for i from #stdPts to r-1 list matrix{for i to n list random(o.IntRange)};
    return join(stdPts,randPts);
)

vanishingIdeal = method()
vanishingIdeal Matrix := pt -> minors(2, pt || basis(1,R))
vanishingIdeal List := pts -> intersect(apply(pts, vanishingIdeal))

evalMatrix = (pts, polys) -> matrix table(pts, polys, (p,f) -> sub(f,p))

degreeGenerators = method(Options => true)
degreeGenerators (List,ZZ) := {IntersectMaxIdeals => true} >> o -> (pts,d) -> (
    if o.IntersectMaxIdeals then
        return super basis(d, vanishingIdeal(pts));
    mons := basis(d,R);
    if #pts == 0 then return mons;
    A := evalMatrix(pts, flatten entries mons);
    return mons * promote(gens ker A,R);
)

IGenD = method()
IGenD (List,ZZ) := (pts,d) -> ideal(degreeGenerators(pts,d))

expectedNumGens = method()
expectedNumGens ZZ := r -> (d := dr(r); {(d, (binomial(n+d,n)-r)), (d+1, max(binomial(n+d+1,n)-r-(n+1)*(binomial(n+d,n)-r), 0))})

isInteresting = method()
isInteresting ZZ := ds -> ds_0_1 > n and ds_1_1 > 0