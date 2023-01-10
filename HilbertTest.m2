restart

n = 2
R = QQ[x_0..x_n]

randomQPoints = method(Options => true)
randomQPoints := {UseStdPoints => true, IntRange => (-42,42)} >> o -> r -> (
    if o.UseStdPoints then (
        stdPts = apply(entries gens QQ^(n+1), e -> matrix{e});
        stdPts = append(stdPts, sum(stdPts));
        if r <= n+2 then return take(stdPts,r);
    ) else stdPts = {};
    randPts = for i from #stdPts to r-1 list matrix{for i to n list random(o.IntRange)};
    return join(stdPts,randPts);
)

evalMatrix = method()
evalMatrix := (pts, polys) -> matrix table(pts, polys, (p,f) -> sub(f,p))

degreeGenerators = method()
degreeGenerators := (pts,d) -> (
    mons := basis(d,R);
    if #pts == 0 then return mons;
    A := evalMatrix(pts, flatten entries mons);
    return mons * promote(gens ker A,R)
)

IgenD = method()
IgenD := (pts,d) -> ideal(flatten entries degreeGenerators(pts,d))

idealFromPoint = method()
idealFromPoint := p -> minors(2, p || basis(1,R))

Jall = method()
Jall := pts -> intersect(apply(pts, idealFromPoint))

d = 5
s_d = binomial(n+d,n)
r = 18
Z = randomQPoints(r)

Mred = R^1/Jall(Z)
reduceHilbert hilbertSeries(Mred)
for i to 10 list hilbertFunction(i, Mred)
res Mred

M = R^1/IgenD(Z,d)
reduceHilbert hilbertSeries(M)
for i to 10 list hilbertFunction(i, M)
res M