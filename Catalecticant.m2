genericPolynomial = (n,d) -> (
    indi := (flatten @@ exponents) \ flatten entries basis(d,QQ[x_1..x_n]);
    R := QQ[indi / (I -> c_I)][x_0..x_(n-1)];
    return sum(indi, I -> c_I*product(gens R, I, (v,i)->v^i));
)

F = genericPolynomial(3,3)

Cat = method()
Cat (ZZ,RingElement) := (d,F) -> (
    dualmons := basis(d,ring F);
    mons := flatten entries basis((degree F)#0 - d, ring F);
    partials = flatten entries contract(dualmons, F);
    return sub(fold( apply(partials,  DF -> last coefficients(DF, Monomials=>mons)), (X,Y) -> X|Y), coefficientRing ring F);
)

I = minors(2, Cat(1, genericPolynomial(3,3)))

leadTerm I