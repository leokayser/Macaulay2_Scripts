restart

n = 3
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

degreeGenerators2 = method()
degreeGenerators2 := (pts,d) -> (
    mons := basis(d,R);
    if #pts == 0 then return mons;
    A := evalMatrix(pts, flatten entries mons);
    return mons * promote(gens ker A,R)
)

idealFromPoint = method()
idealFromPoint := p -> minors(2, p || basis(1,R))

Jall = method()
Jall := pts -> intersect(apply(pts, idealFromPoint))

degreeGenerators1 = method()
degreeGenerators1 := (pts,d) -> super basis(d, Jall(pts))

IgenD = method()
IgenD := (pts,d) -> ideal(degreeGenerators1(pts,d))


d = 5
s_d = binomial(n+d,n)
r = s_d -n - 1
Z = randomQPoints(r)
I = IgenD(Z,d);

M = R^1/I;
for i to 15 list hilbertFunction(i, M)

--elapsedTime basis(d,Jall(Z));
--elapsedTime degreeGenerators(Z,d);

-- Mred = R^1/Jall(Z);
-- reduceHilbert hilbertSeries(Mred)
-- for i to 10 list hilbertFunction(i, Mred)
-- res Mred

elapsedTime M = R^1/IgenD(Z,d);
reduceHilbert hilbertSeries(M)
for i to 20 list hilbertFunction(i, M)
res M