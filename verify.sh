#! /bin/bash
cmda="load \"TestSeq.m2\";testSpecific($1,$2,IdealMethod=>$3);"
echo $cmd
M2 <<< $cmda &
disown