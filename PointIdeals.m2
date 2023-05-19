n = 3
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

hS = d -> binomial(n+d,n)

firstHFeq = method()
firstHFeq ZZ := r -> (d := 1; while binomial(n+d,n) < r do (d = d+1); return d)

expectedNumGens = method()
expectedNumGens ZZ := r -> (d := firstHFeq(r); {(d, hS(d)-r), (d+1, max(hS(d+1)-r-(n+1)*(hS(d)-r), 0))})

isInteresting = method()
isInteresting ZZ := r -> (g := expectedNumGens(r); g_0_1 > n and g_1_1 > 0)

expectedHilbFun = method()
expectedHilbFun (ZZ, List, ZZ) := (r,g,i) -> (
    d := g_0_0;
    if i < d then return hS(i);
    return max(hS(i)-g_0_1*hS(i-d), r)
)
expectedHilbFun (ZZ,ZZ) := (r,i) -> expectedHilbFun(r, expectedNumGens(r), i)


