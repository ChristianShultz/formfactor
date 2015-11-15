#!/bin/tcsh -x

# copy this guy down 
set ini = check_rots.@@@.ini.xml

set badlist = bad-Q2s.list

if ( ! -f $ini ) then 
  echo "missing ini: $ini"
  exit 1
endif

if ( -f $badlist) then 
  echo "removing bad list"
  /bin/rm $badlist
endif

# we are running the single Q2 variant of radmat
set exen = radmat_util
set exe = `which $exen`

if ( ! -e $exe ) then 
  echo "$exen : $exe is not an executable"
  exit 1
endif

# find all of the Q2 dirs
foreach dir ( `find -maxdepth 2 -type d -ls | awk '{print $11}' | grep lefty | cut -c 3-` )
  
  # what are we doing
  echo "working on $dir"
    
  pushd $dir
  pushd llsq

  set sdb = state_database.rad
  set ros = row_index_to_continuum_elem.txt
  set xml = checkrot.xml

  echo " cat ../../../$ini > $xml " | $SHELL -x

  $exe rot_llsq $xml 1 1 

  popd 
  popd

end

