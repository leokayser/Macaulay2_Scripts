load "Chopnumbers.m2"

n = 2
p = 1009
savePoints = true

testAll = method()
testAll (ZZ,ZZ) := (dmin,dmax) -> (

    S = ZZ/p[x_0..x_n];
    for d from dmin to dmax do (
        print("-------- d = " | d | " --------");

        (rmin, rmax) = interestingRange(n,d);
        print("--- r = " | rmin | " ... " | rmax | " ---");

        for r from rmin to rmax do (
            if r == rmax then (
                print("(d=" | d | ",rmax=" | rmax | ") is already proved.");
                continue;
            );
            elapsedTime (pts := testCase(S,r));
            if savePoints then (
                file = "pts/P" | n | "_r" | r | ".txt";
                file << toString (p,pts) << close;
            );
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