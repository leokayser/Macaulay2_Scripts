restart

load "HilbertTest.m2"

d = 3
a = 0..(d-1)
pts = apply(toList((0,0)..(n,n)), t -> matrix{{1,a_(t_0),a_(t_1)}})

degreeGenerators1(pts,d)