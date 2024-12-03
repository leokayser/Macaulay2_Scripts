R = QQ[x,y,z];
A = {x,y,z,x+y,x+z,y+z}

needsPackage "HyperplaneArrangements";
M = der arrangement A

R = QQ[x,y,z,s1,s2,s3,s4,s5,Degrees => {{1,0},{1,0},{1,0},{0,1},{0,1},{0,1},{0,1},{0,1}}];

s6 = -(s1+s2+s3+s4+s5);
I = ideal(s1+s2+s3+s4+s5+s6,
     x*y*s3+x*z*s3+y*z*s3+z^2*s3+y*z*s5+z^2*s5+x*z*s6+z^2*s6,
     x*y*s2+y^2*s2+x*z*s2+y*z*s2+y^2*s4+y*z*s4+x*y*s6+y^2*s6,
     x^2*s2-y^2*s2+x^2*s3-z^2*s3+x*y*s4-y^2*s4+x*z*s5-z^2*s5+x^2*s6-y^2*s6+y*z*s6-z^2*s6);
multidegree I

f = first first entries gens eliminate({z,y},I+ideal(y-1)); d1 = discriminant(f,x);
f = first first entries gens eliminate({y,x},I+ideal(x-1)); d2 = discriminant(f,z);
f = first first entries gens eliminate({x,z},I+ideal(z-1)); d3 = discriminant(f,y);

dis = gcd(d1,d2,d2);
toString factor(dis)
# terms dis
