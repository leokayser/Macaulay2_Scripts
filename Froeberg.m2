test = (p,ds) -> (
    n := #ds - 1;
    kk := if p == 0 then QQ else ZZ/p;
    S := kk[x_1..x_n];
    fs := apply(gens S, ds_{0..(n-1)}, (z,i)->z^i) | {(sum gens S)^(ds#n)};
    I := ideal(fs);
    return S/I;
)

ratHS = M -> (
    hs := reduceHilbert hilbertSeries M;
    num := (value numerator hs);
    den := (value denominator hs);
    return if num % den == 0 then num // den else num / den;
)

expectedHS = (n,ds) -> (
    A := ZZ[x];
    use degreesRing A;
    num := product(ds_{0..(min(n,#ds)-1)}, d -> 1-T^d);
    den := (1-T)^n;
    if #ds < n then return num/den;
    p := num // den * product(ds_{n..(#ds-1)}, d -> 1-T^d);
    return fceil p;
)

expectest = ds -> expectedHS(#ds - 1, ds);

coeffList = p -> (
    T := (ring p)_0;
    return for i to (degree p)_0 list coefficient(T^i, p);
)

fceil = p -> (
    L := coeffList p;
    T := (ring p)_0;
    e := #L - 1;
    for i to e do (
        if L_i <= 0 then (e = i-1; break;);
    );
    L' := L_{0..e};
    return sum(0..e, L', (i,c) -> c*T^i);
)

equals = method();
equals (Divide,Divide) := (s1,s2) -> (
    s2' := sub(s2, ring numerator s1);
    den1 := denominator s1;
    num1 := numerator s1;
    den2 := denominator s2';
    num2 := numerator s2';
    return value( den1*num2 ) == value( den2*num1 );
)