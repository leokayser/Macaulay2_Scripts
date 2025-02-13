restart
kk = ZZ/101; KK = frac kk[t]; S = KK[x_0..x_3]
G = (x_0^2 + t*x_1^2)^2 + t^3*x_2^4
P = G^3 + G*t^7*x_3^8

d = 12

veronese = (d,S) -> (
    mons = flatten entries basis(d,S);
    R := (coefficientRing S)[T_1..T_#mons];
    return map(S,R,mons);
)
phi = veronese(12,S);
ker phi
