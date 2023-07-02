#! /bin/bash
cmda="load \"TestSeq.m2\";n=$1;testAll($2,$3,IdealMethod=>f1);"
cmdb="load \"TestSeq.m2\";n=$1;testAll($2,$3,IdealMethod=>f2);"
echo $cmd
M2 <<< $cmda &
disown
M2 <<< $cmdb &
disown