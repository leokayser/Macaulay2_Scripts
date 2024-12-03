--load "Chopnumbers.m2"

mu1 = (n,d,r) -> (
    gend := max(binomial(n-1+d,n-1)-r,0);
    S1Id := gend*n;
    Id1 := max(binomial(n+d,n-1)-r,0);
    return min(S1Id,Id1);
)

edim = (n,d,r) -> max(binomial(n-1+d,n-1)-r,0)*n

delta = (n,d,r) -> edim(n,d,r) - mu1(n,d,r)

delta' = (n,d,r) -> delta(n,d,r) - delta(n,d,r+1)

rmax = (n,d) -> binomial(n-1+d,n-1)   -- Here I(p1..pr)_d = 0 and all generators in d+1, none in d
rmin = (n,d) -> binomial(n-1+d-1,n-1) -- Here I(p1..pr)_(d-1) = 0 and all generators in d, hence mult. is surjective

(n,d) = (5+1,3)
netList transpose ({{"r","d","d'"}} | for i from 0 to rmax(n,d)-rmin(n,d) list (r = rmax(n,d)-i; {r,delta(n,d,r),delta'(n,d,r)}))