load "Chopnumbers.m2"

n = 2
p = 1009
savePts = true

padInt = (n,l) -> (
    str := toString(n);
    while #str < l do str = "0" | str;
    return str;
)

testCase = method(Options => true)
testCase (PolynomialRing,ZZ) := {SavePoints => false} >> o -> (S,r) -> (
    kk := coefficientRing S;
    n := (numgens S)-1;
    d := regH(n,r);

    identStr := concatenate("(n=", toString(n), ",d=", toString(d) , ",r=", toString(r), ")");
    print(identStr | " -> expected gap = " | expectedGapSize(n,r));
    
    while true do (
        randPts := randomPoints_S(r);
        I := vanishIdeal_S(randPts);
        if verifyConj(I,r) then (
            print(identStr | " valid!");
            if o.SavePoints then (
                file = "pts/P" | n | "_d" | padInt(d,2) | "_r" | padInt(r,3) | ".txt";
                file << toString (n,p) << endl;
                for pt in randPts do (file << toString(pt) << endl);
                file << close;
            );
            return randPts;
        );
        print(identStr | " invalid, trying again");
    );
)


testAll = method()
testAll (ZZ,ZZ) := (dmin,dmax) -> (

    S = ZZ/p[x_0..x_n];
    for d from dmin to dmax do (
        print("-------- d = " | d | " --------");

        (rmin, rmax) = interestingRange(n,d);
        print("--- r = " | rmin | " ... " | rmax | " ---");

        for r from rmin to rmax do (
            if r == rmax then (
                print("(d=" | d | ",rmax=" | rmax | ") is already proven.");
                continue;
            );
            elapsedTime (pts := testCase(S,r,SavePoints=>savePts));
        );
    );
)

-- testAll = method()
-- testAll (ZZ,ZZ) := (dmin,dmax) -> (
--     tasks = {};
--     for d from dmin to dmax do (
--         (rmin, rmax) = interestingRange(n,d);
--         tasks = join(tasks, for r from rmin to rmax list createTask((d,r) -> testCase(n,d,r), (d,r)));
--     );
--     numTasks := #tasks;
--     print("ToDo: " | numTasks | " cases");
--     apply(tasks, schedule);
--     while #tasks > 0 do (
--             scan(tasks, t -> if isReady(t) then tasks = delete(t,tasks));
--             sleep 1;
--         );
--     print("All " | numTasks | " cases done!");
-- )