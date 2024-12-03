ecompi = method();
ecompi (RingElement,RingElement,Ideal) := (F,t,I) -> (
    S := ring t;
    e := (degree t)_0;
    assert(number(flatten entries mingens I, f -> (degree f) != {e}) == 0);
    Fperp := inverseSystem F;
    d := length module (S/((Fperp : I) + ideal(t)));
    assert(d < infinity);
    return ceiling (d/e);
)
ecomp = method();
ecomp (RingElement,RingElement) := (F,t) -> ecompi(F,t,ideal t);



--KK = toField(QQ[i]/ideal(i^2+1));
makeQuad = (k,l) -> (
    use QQ[x_1..x_k,y_1..y_l];
    F := sum(1..(k//2),i->x_(2*i) * x_(2*i-1));
    G := sum(1..(l//2),i->y_(2*i) * y_(2*i-1));
    if k % 2 == 1 then F = F + x_k^2;
    if l % 2 == 1 then G = G + y_l^2;
    return (F,G);
)

n = 4
(F,G) = makeQuad (4,1)
I = ideal(x_1..x_n,y_1..y_1)
t = x_1
ecompi(F*G, t, I)

Fperp = inverseSystem(F*G)
flatten entries mingens (Fperp : I)