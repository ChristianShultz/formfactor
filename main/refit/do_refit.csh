#!/bin/tcsh 

# copy this guy down 
set ini = refit.@@@.ini.xml

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
foreach dir ( `ls -al | grep ^d  | grep -v refit | grep -v FF | awk '{print $9}' | grep Q2` )

    set refit = $dir.refit
    mkdir $refit 

    # push refit
    pushd $refit 

    foreach group( `ls ../$dir | grep -v continuum | grep -v correlator` )

      echo "working on $dir/$group \n"
      set working = ../$dir/$group
      set llsq = $working/llsq

      # break if no llsq dir 
      if ( ! -d $llsq ) then 
        echo "$llsq is missing \n"
        continue 
      endif 

      # push group 
      mkdir $group
      pushd $group     

      # find files and set soft links
      set llsqsdb = ../$llsq/state_database.rad
      set sdb = state_database.rad 
      set ffsdb = ../$working/ff_database.rad
      set ff = ff_database.rad
      set llsqelems = ../$llsq/row_index_to_continuum_elem.txt
      set elems = row_to_cont.txt

      if ( ! -f $llsqsdb ) then 
        echo "$llsqsdb is missing "
        continue
      endif

      if ( ! -f $llsqelems ) then 
        echo "$llsqelems is missing"
        continue 
      endif 

      if ( ! -f $ffsdb ) then 
        echo "$ffsdb is missing"
        continue 
      endif 

      ln -s $llsqsdb $sdb 
      ln -s $llsqelems $elems 
      ln -s $ffsdb $ff 

      # find the sequence of elems 
      set lats = ` cat $elems | awk '{ print $1 }' | xargs `
      
      # make an ini via sed 
      set xml = ini.xml
      
      # sed to sub in on the ini
      echo " cat ../../$ini | sed 's/@@@@/$lats/g' > $xml " | $SHELL -x

      $exe refit_ffs $xml  

      # pop group 
      popd
      
    end
    
    #pop refit 
    popd 

end

