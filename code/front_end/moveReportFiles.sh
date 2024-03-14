#!/bin/bash

#!/bin/bash

protname="$1"

alreadyrun="../proteinsparams/listoffiles.txt"
datafolder="../protein"

found=`egrep ^$protname, $alreadyrun`
if [ ! -z "$found" ]; then
        #paramfile=`egrep ^$protname $alreadyrun | awk -F',' '{print $2}'`
        protfolder=`egrep ^$protname, $alreadyrun | awk -F',' '{print $2}'`
        
        mv "$protfolder"* $datafolder"/"$protfolder"/progress_msgs/."

else
	protfolder=$protname"_allparam"	
        mv "$protname"* $datafolder"/"$protfolder"/progress_msgs/."
fi

        
