load "Chopnumbers.m2"
needsPackage("Points")

f1 = (ptlst,S) -> (
    pts1 := matrix copy ptlst;
    R1 := coefficientRing(S)[x_0..x_(numgens S - 1)];
    use R1;
    mat := sub(pts1,R1);
    return points(mat);
)
f2 = (ptlst,S) -> (
    pts2 := matrix copy ptlst;
    R2 := coefficientRing(S)[x_0..x_(numgens S - 1)];
    use R2;
    (inG,G) := projectivePoints(pts2,R2,VerifyPoints=>false);
    return ideal G;
)

-- Uses parallel tasks to find the ideal of points
tryBothMethods = method()
tryBothMethods (Matrix,Ring) := (pts,S) -> ( 
    allowableThreads = 3;
    t1 := schedule(f1,(entries pts,S));
    t2 := schedule(f2,(entries pts,S));
    while true do (
        nanosleep 1000000000;
        if isReady(t1) then (
            I1 := taskResult(t1);
            if instance(I1,Nothing) then(logProgress("points error!"); continue;);
            logProgress("points won!");
            try (if isReady(t2) then taskResult(t2) else (try cancelTask(t2); print("t2 terminated")));
            return I1;
        );
        if isReady(t2) then (
            I2 := taskResult(t2);
            if instance(I2,Nothing) then(logProgress("projectivePoints error!"); continue;);
            logProgress("projectivePoints won!");
            try (if isReady(t1) then taskResult(t1) else (try cancelTask(t1); print("t1 terminated")));
            return I2;
        );
    );
)


p = 1009


testCase = method()
testCase (ZZ,ZZ) := (n,r) -> (
    S := ZZ/p[vars(0..n)];
    d := regH(n,r);

    identStr := concatenate("(n=", toString(n), ",d=", toString(d) , ",r=", toString(r), ")");
    print(identStr | " -> expected gap = " | expectedGapSize(n,r));
    
    while true do (
        randPts := random(ZZ^(n+1),ZZ^r, Height=>p);
        elapsedTime (I := ideal((projectivePoints(randPts,S))_1));
        print("Ideal calculated.");
        if verifyESC(I,r) then (
            print(identStr | " valid!");
            return randPts;
        );
        print(identStr | " invalid, trying again");
    );
)



padInt = (n,l) -> (
    str := toString(n);
    while #str < l do str = "0" | str;
    return str;
)
filename = (dir,n,d,r) -> concatenate(dir,"/P",toString(n),"_d",padInt(d,2),"_r",padInt(r,3),".txt");
betterfilename = (n,r) -> concatenate("P",padInt(n,2),"/d",padInt(regH(n,r),2),"_r",padInt(r,4),".txt");

savePoints = method()
savePoints (Matrix,String) := (pts,file) -> (
    ptsList := entries transpose pts;
    file << toString (n,p) << endl;
    for pt in ptsList do (file << toString(pt) << endl);
    file << close;
)
betterSavePoints = method()
betterSavePoints (Matrix,Ideal,String) := (pts,I,file) -> (
    ptsList := entries transpose pts;
    file << toExternalString ring I << endl;
    file << toExternalString (#ptsList) << endl;
    for pt in ptsList do (file << toExternalString(pt) << endl);
    d := regH(#(ptsList_0) - 1,#ptsList);
    gs := select(first entries gens I, f -> (degree f)_0 == d);
    for g in gs do file << toExternalString(g) << endl;
    file << close;
)

timestamp = () -> (
    t := currentTime() + 2*(3600); -- Hard coded time diff.
    return concatenate(
        toString(t//(60*60*24)),"-",
        padInt((t//(60*60))%24,2),":",
        padInt((t//60)%60,2),":",
        padInt(t%60,2),"> ");
)
logProgress = msg -> (
    t := currentTime();
    f := openOutAppend("log.txt");
    print(timestamp() | msg);
    f << timestamp() << msg << endl << close;
)

loadPoints = method()
loadPoints (String) := file -> (
    f := lines get file;
    (n,p) := value f#0;
    return transpose matrix for i from 1 to #f - 1 list value f#i;
)
betterLoadPoints = method()
betterLoadPoints (String) := file -> (
    f := lines get file;
    S := value f#0;
    r := value f#1;
    pts := matrix transpose for i from 2 to r+1 list value f#i;
    I := ideal for i from r+2 to #f - 1 list value f#i;
    return (pts,I);
)

transformAll = method()
transformAll (ZZ) := n -> (
    prefix := "P" | padInt(n,2);
    if not fileExists prefix then mkdir prefix;
    S := ZZ/p[x_0..x_n];
    prevCalcs := sort(readDirectory("points"));
    logProgress("Starting transform "|toString(n));
    for file in prevCalcs do (
        if not match("P"|toString(n), file) then continue;
        --match("r.*",file);
        --r := value(substring((lastMatch_0_0 + 1,lastMatch_0_0 - 3), file));
        --print(class(r));
        pts := loadPoints("points/"|file);
        r := numgens source pts;
        newfile := betterfilename(n,r);
        if fileExists(newfile) then continue;
        logProgress("Loaded file "|file);
        I := tryBothMethods(pts,S);
        logProgress(newfile | " ideal calculated");
        assert verifyESC(I,r);
        logProgress(newfile | " ESC verified");
        betterSavePoints(pts,I,newfile);
    );
    logProgress(">>>>>>>>>>        "|toString(n)|" DONE!        <<<<<<<<<<");
    exit;
)

testAll = method()
testAll (ZZ,ZZ) := (dmin,dmax) -> (

    igc := if n <= 4 then true else false;

    for d from dmin to dmax do (
        print("-------- d = " | d | " --------");
        (rmin, rmax) := interestingRange(n,d, Igc=>igc);
        print("--- r = " | rmin | " ... " | rmax | " ---");

        for r from rmin to rmax do (
            if r == rmax then (
                print("(d=" | d | ",rmax=" | rmax | ") is already proven.");
                continue;
            );
            file := filename("points",n,d,r);
            if fileExists(file) then (
                print(concatenate("\"",file,"\" exists."));
                continue;
            );
            elapsedTime (pts := testCase(n,r));
            savePoints(pts,file);
        );
    );
)
