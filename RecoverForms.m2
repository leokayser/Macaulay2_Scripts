HFtable = method(Dispatch => Thing, Options => {Range => null})
HFtable Sequence := o -> Is -> (
    idealList = toList Is;
    d := if o.Range === null then max apply(idealList, I -> regularity I) else o.Range;
    netList transpose for i from 0 to d list {i} | apply(idealList, I -> hilbertFunction({i},I))
)
HFtable Ideal := o -> I -> (
    d := if o.Range === null then regularity I else o.Range;
    netList transpose for i from 0 to d list {i, hilbertFunction({i},I)}
)

degStr := f -> if f == 0 then " " else toString sum degree f;

degTable = method()
degTable Matrix := M -> netList applyTable(entries M, degStr)
degTable ChainComplex := C -> apply(1..(length C), i -> degTable C.dd_i)

chop = method()
chop (Ideal,ZZ) := (I,d) -> ideal select(flatten entries gens I, f -> (degree f)#0 == d)
chop Ideal := I -> (
    d := min apply(flatten entries gens I, f -> (degree f)#0);
    chop(I,d)
)


end
--------
restart
load "RecoverForms.m2"

S = QQ[x_0..x_2]
d = 5; t = (d+1)/2;
F = random(d, S, Height=>2); I = inverseSystem(F);
minimalBetti I
HFtable(I, chop I)
degTable res I


q = for i from 1 to 6 list random(2,S, Height=>2)
auxRing = QQ[c_0..c_9]
P = sub(genericSkewMatrix(auxRing,5), matrix{{x_0,x_1,0,0,q_0,q_1,q_2,q_3,q_4,q_5}});

assert(P + transpose P == 0); degTable P
I = pfaffians(4,P)
HFtable(I, chop I)