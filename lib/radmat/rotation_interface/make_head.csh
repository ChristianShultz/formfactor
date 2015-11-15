#!/bin/tcsh 

set pth = `pwd`
set f = `basename $pth` 
set ff = $f.h 

set hg = ` echo $f | tr "[a-z]" "[A-Z]"`
set guard = ${hg}_H

echo "#ifndef $guard" > $ff
echo "#define $guard" >> $ff

ls *.h | grep -v $ff |  awk '{ print "#include \"" $1 "\"" }' >>  $ff 

echo "#endif" >> $ff

