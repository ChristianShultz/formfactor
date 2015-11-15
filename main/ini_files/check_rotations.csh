#!/bin/tcsh -x 

set out =  check_rots.out 
nohup ./check_rots_llsq.csh >& $out & 
