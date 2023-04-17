#!/usr/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/star/u/rcorliss/install/
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
 
for (( n=0; n<30; n++ ));
do 
    echo Processing $1 $2 $3 with $4 divisions.  thisjob is $5 ;
    root -b -q ./generate_rossegger.C\($1,$2,$3,$4,$5\)
done

echo all done
