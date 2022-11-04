#!/bin/bash
#cd CORE
#./RACIPE *.topo -num_paras 7000 &
for i in {1..6}
	do 
		mkdir OEid$i
		cd OEid$i
		scp ../../../RACIPE_dummy/* .
		scp ../../../core.topo .
		./RACIPE core.topo -num_paras 7000 -OEID $i -OEFD 30 
		cd ..
	done
