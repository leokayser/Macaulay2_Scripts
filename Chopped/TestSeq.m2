load "Chopnumbers.m2"
needsPackage("Points")

f1 = (pts,S) -> points(sub(pts,S));
f2 = (pts,S) -> ideal last projectivePoints(pts,S,VerifyPoints=>false);

p = 1009

testCase = method(Options => {IdealMethod => f1})
testCase (ZZ,ZZ) := o -> (n,r) -> (
    S := ZZ/p[x_0..x_n];
    d := regH(n,r);

    identStr := concatenate("(n=", toString(n), ",d=", toString(d) , ",r=", toString(r), ")");
    print(identStr | " -> expected gap = " | expectedGapSize(n,r));
    
    newfile := betterfilename(n,r);

    while true do (
        randPts := random(ZZ^(n+1),ZZ^r, Height=>p);
        elapsedTime (I := o.IdealMethod(randPts,S));
        if fileExists(newfile) then return null;
        logProgress(betterfilename(n,r)|" Ideal calculated by "|toString(o.IdealMethod));
        if verifyESC(I,r) then (
            if fileExists(newfile) then(logProgress(betterfilename(n,r)|" "|toString(o.IdealMethod)|" skip."); return null;);
            print(identStr | " valid!");
            return (randPts,I);
        );
        if fileExists(newfile) then return null;
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

--savePoints = method()
--savePoints (Matrix,String) := (pts,file) -> (
--    ptsList := entries transpose pts;
--    file << toString (n,p) << endl;
--    for pt in ptsList do (file << toString(pt) << endl);
--    file << close;
--)
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

transformAll = method(Options => {IdealMethod => f1})
transformAll (ZZ) := o -> n -> (
    prefix := "P" | padInt(n,2);
    if not fileExists prefix then mkdir prefix;
    S := ZZ/p[x_0..x_n];
    prevCalcs := sort(readDirectory("points"));
    logProgress("Starting transform "|toString(n));
    for file in prevCalcs do (
        if not match("P"|toString(n), file) then continue;
        match("r.*",file);
        m := lastMatch_0_0;
        r := value substring((m + 1,#file-m-5), file);
        newfile := betterfilename(n,r);
        if fileExists(newfile) and not fileExists(newfile | toString(o.IdealMethod)) then continue;
        pts := loadPoints("points/"|file);
        logProgress("Loaded file "|file);
        I := o.IdealMethod(pts,S);

        if fileExists(newfile) and not fileExists(newfile | toString(o.IdealMethod)) then
            (logProgress("A: "|toString(o.IdealMethod)|" skipping ahead."); continue;);
        if toString(o.IdealMethod) == "f1" and fileExists(newfile | "f2") then
            (logProgress("B: f1 skipping ahead."); continue;);
        if toString(o.IdealMethod) == "f2" and fileExists(newfile | "f1") then
            (logProgress("B: f2 skipping ahead."); continue;);

        (newfile | toString(o.IdealMethod)) << endl << close; -- Signal computation in progress
        logProgress(newfile | " ideal calculated");
        betterSavePoints(pts,I,newfile);
        if verifyESC(I,r) then (
            logProgress(newfile | " ESC verified");
            removeFile(newfile | toString(o.IdealMethod));
        ) else logProgress(newfile | " ESC not verified !!!!!!!!!!!!!!!");
    );
    logProgress(">>>>>>>>>>        "|toString(n)|" DONE!        <<<<<<<<<<");
    quit;
)

testSpecific = method(Options => {IdealMethod => f1})
testSpecific (ZZ,ZZ) := o -> (n,r) -> (
    d := regH(n,r);
    file := betterfilename(n,r);
    if fileExists(file) then (
        logProgress(concatenate("\"",file,"\" exists."));
        quit;
        return;
    );
    result := testCase(n,r, IdealMethod=>o.IdealMethod);
    if instance(result,Nothing) then quit;
    logProgress(file | " done.");
    (pts,I) := result;
    betterSavePoints(pts,I,file);
    quit;
)

testAll = method(Options => {IdealMethod => f1})
testAll (ZZ,ZZ) := o -> (dmin,dmax) -> (
    prefix := "P" | padInt(n,2);
    if not fileExists prefix then mkdir prefix;
    igc := if n <= 4 then true else false;

    for d from dmin to dmax do (
        --print("-------- d = " | d | " --------");
        (rmin, rmax) := interestingRange(n,d, Igc=>igc);
        --print("--- r = " | rmin | " ... " | rmax | " ---");

        for r from rmin to rmax do (
            if r == rmax then (
                --print("(d=" | d | ",rmax=" | rmax | ") is already proven.");
                continue;
            );
            file := betterfilename(n,r);
            if fileExists(file) then (
                --print(concatenate("\"",file,"\" exists."));
                continue;
            );
            result := testCase(n,r, IdealMethod=>o.IdealMethod);
            if instance(result,Nothing) then continue;
            logProgress(file | " done.");
            (pts,I) := result;
            betterSavePoints(pts,I,file);
        );
    );
    logProgress(">>>>>  P^" | toString(n)| " all degrees [" | toString(dmin) | "," | toString(dmax) | "] verified!  <<<<<");
    quit;
)
