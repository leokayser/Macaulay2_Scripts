r = 3
R = QQ[a_0..a_r]


T = k -> (
    Toep := mutableMatrix(R, k, r+k);
    aa := (for i from 0 to r list a_(r-i));
    for i from 0 to k-1 do for j from 0 to r do Toep_(i,i+j) = aa#j;
    return matrix Toep;
)

A = T(3)
I = minors(3,A)

-- T(m+1) has size (m+1) x (r+m+1), rank m locus
-- expected codimension is (m+1-m)(r+m+1-m) = r+1
-- V(minors) = empty, has dimension -1, codim_P^r = r+1

-- S = K[a0..ar], dim S_d = binomial(d+r, r)
-- I = <S_{m+1}> if and only if deg M_m(T) = dim S/I = sum dim S_i = binomial(r+1+m, r+1)

-- Apply 3264 Theorem 12.5:
-- e = m+1
-- f = r+m+1
-- binomial(r+m+1, r+1)

dim W = r+1
P(W) -> P(Sym^{m+1}W)
W = S^r(V^*), V^* = K[x,y]_1
P(S_r) -> P(Sym^{m+1}S_r)

P(S_r) -> Gr(m+1, S_{r+m}) -> P(A^{m+1} S_{r+m})