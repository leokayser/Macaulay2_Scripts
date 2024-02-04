needsPackage "SubalgebraBases";

homogeneousBinary = n -> (
    use QQ[L,R,O_1..O_n,I_1..I_n,H_0..H_n]; --L,R states; O_i,I_i cells; H_i head
    
    rightTrans0 := for i from 1 to n   list -R*O_i*H_i + L*I_i*H_(i-1);
    rightTrans1 := for i from 1 to n-1 list -R*I_i*H_i + R*O_i*H_(i+1);
    leftTrans   := for i from 1 to n   list -L*O_i*H_i + L*O_i*H_(i-1);
    
    allTrans := rightTrans0 | rightTrans1
              | leftTrans | {-L*H_0 + R*H_1}; --Turn around on the left

    padVars  := toList(O_1..O_n) | toList(I_1..I_n);

    f := L*H_0*product(toList(O_1..O_n))
       - R*H_n*product(toList(O_1..O_(n-1)))*I_n;

    return (subring(allTrans | padVars), f);
)

homogeneousBinary = n -> (
    use QQ[L,R,O_1..O_n,I_1..I_n,H_0..H_n]; --L,R states; O_i,I_i cells; H_i head
    
    rightTrans0 := for i from 1 to n   list -R*O_i*H_i + L*I_i*H_(i-1);
    rightTrans1 := for i from 1 to n-1 list -R*I_i*H_i + R*O_i*H_(i+1);
    leftTrans   := for i from 1 to n   list -L*O_i*H_i + L*O_i*H_(i-1);
    
    allTrans := rightTrans0 | rightTrans1
              | leftTrans | {-L*H_0 + R*H_1}; --Turn around on the left

    padVars  := toList(O_1..O_n) | toList(I_1..I_n);

    f := L*H_0*product(toList(O_1..O_n))
       - R*H_n*product(toList(O_1..O_(n-1)))*I_n;

    return (subring(allTrans | padVars), f);
)

(A,f) = homogeneousBinary(3)
groebnerMembershipTest(f,A)
f % A
f // A