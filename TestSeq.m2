load "Chopnumbers.m2"
loadPackage "Points"

n = 8
p = 1009


testCase = method()
testCase (ZZ,ZZ) := (n,r) -> (
    S := ZZ/p[vars(0..n)];
    d := regH(n,r);

    identStr := concatenate("(n=", toString(n), ",d=", toString(d) , ",r=", toString(r), ")");
    print(identStr | " -> expected gap = " | expectedGapSize(n,r));
    
    while true do (
        --randPts := randomPointsMat(S,r);
        --elapsedTime (I := points(randPts));
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

savePoints = method()
savePoints (Matrix,String) := (pts,file) -> (
    ptsList = entries transpose pts;
    file << toString (n,p) << endl;
    for pt in ptsList do (file << toString(pt) << endl);
    file << close;
)

testAll = method()
testAll (ZZ,ZZ) := (dmin,dmax) -> (

    igc = if n <= 4 then true else false;

    for d from dmin to dmax do (
        print("-------- d = " | d | " --------");
        (rmin, rmax) = interestingRange(n,d, Igc=>igc);
        print("--- r = " | rmin | " ... " | rmax | " ---");

        for r from rmin to rmax do (
            if r == rmax then (
                print("(d=" | d | ",rmax=" | rmax | ") is already proven.");
                continue;
            );
            file = filename("points",n,d,r);
            if fileExists(file) then (
                print(concatenate("\"",file,"\" exists."));
                continue;
            );
            elapsedTime (pts := testCase(n,r));
            savePoints(pts,file);
        );
    );
)