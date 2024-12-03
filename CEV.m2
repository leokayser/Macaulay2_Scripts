restart
S = QQ[x_1..x_3]
d = 3; t = (d+1)/2;
F = random(d,S, Height=>2); I = inverseSystem(F);
minimalBetti I
A = (res I).dd_2;
A' = fold((M,N)->M||N, for i from 0 to 4 list A^{4-i} )

deg = f -> if f == 0 then " " else sum degree f
netList applyTable(entries A', deg)

d = 3
F = random(d,S, Height=>2); I = inverseSystem(F);
lowdeg = f -> degree f == {(d+1)//2}
J = ideal select(flatten entries gens I, lowdeg);
--minimalBetti J
netList transpose for i from 0 to d+2 list {i,hilbertFunction({i},I),hilbertFunction({i},J)}

-- 0 <- S/I <- S <- S[-3]^16 <- S[-4]^30 <- S[-5]^16 <- S[-8] <- 0

S = QQ[x_1..x_2]
for d from 2 to 10 do (
    print "d=" | toString(d);
    F = random(d,S, Height=>2); I = inverseSystem(F);
    print minimalBetti I;
)


c = for i from 1 to 6 list random(2,S)
P = matrix {{0, c_0, c_1, c_2, x_0},
        {-c_0, 0, c_3, c_4, x_1},
        {-c_1, -c_3, 0, c_5, 0},
        {-c_2, -c_4, -c_5, 0, 0},
        {-x_0, -x_1, 0, 0, 0}}
I = pfaffians(4,P)

I = inverseSystem (x_0^4+x_1^4+x_2^4)

Lemma: If k >= max {j | beta_1,j != 0}, then <DkF> recovers F
Lemma: If beta_1,d != 0, then <DkF> for k<d never recovers F, in fact, the fiber is at least a PP^beta_1,d
