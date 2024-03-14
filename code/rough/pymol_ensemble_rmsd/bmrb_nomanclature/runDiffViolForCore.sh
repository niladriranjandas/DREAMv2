#!/bin/bash

#2m4k

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2m4k/output_pdb/grp_2m4k.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 2m4k_anchor_pseudo_cluster.pdb $resis > 2m4k_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 2m4k_anchor_pseudo_core.pdb 2m4k.dist.1.upl > 2m4k_viol_diffbound_anchor_core.txt

resis=""

#2kok

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2kok/output_pdb/grp_2kok.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 2kok_anchor_pseudo_cluster.pdb $resis > 2kok_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 2kok_anchor_pseudo_core.pdb 2kok.dist.1.upl > 2kok_viol_diffbound_anchor_core.txt

resis=""

#5x1x

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/5x1x/output_pdb/grp_5x1x.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 5x1x_anchor_pseudo_cluster.pdb $resis > 5x1x_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 5x1x_anchor_pseudo_core.pdb 5x1x.dist.1.upl > 5x1x_viol_diffbound_anchor_core.txt

resis=""

#1pbu

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/1pbu/output_pdb/grp_1pbu.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 1pbu_anchor_pseudo_cluster.pdb $resis > 1pbu_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 1pbu_anchor_pseudo_core.pdb 1pbu.dist.1.upl > 1pbu_viol_diffbound_anchor_core.txt

resis=""

#1xxe

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/1xxe/output_pdb/grp_1xxe.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 1xxe_anchor_pseudo_cluster.pdb $resis > 1xxe_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 1xxe_anchor_pseudo_core.pdb 1xxe.upl > 1xxe_viol_diffbound_anchor_core.txt

resis=""

#2lav

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2lav/output_pdb/grp_2lav.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 2lav_anchor_pseudo_cluster.pdb $resis > 2lav_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 2lav_anchor_pseudo_core.pdb 2lav.upl > 2lav_viol_diffbound_anchor_core.txt

resis=""

#2l7b

resis=`./getResiFromPDB.sh /home/niladri/Documents/Disco_etc_all_in_1/our_algo/archive/2l7b/output_pdb/grp_2l7b.pdb`
if [ ! -z $resis ]; then
	./writePDBforRange.sh 2l7b_anchor_pseudo_cluster.pdb $resis > 2l7b_anchor_pseudo_core.pdb
fi

./divBoundsgetViol.sh 2l7b_anchor_pseudo_core.pdb 2l7b.upl > 2l7b_viol_diffbound_anchor_core.txt

resis=""
