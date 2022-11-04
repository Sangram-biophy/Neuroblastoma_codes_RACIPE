#!/bin/bash
#cd CORE
#./RACIPE *.topo -num_paras 7000 &
for i in {1..6}
	do 
		mkdir DEid$i
		cd DEid$i
		scp ../../../RACIPE_dummy/* .
		scp ../../../core.topo .
		./RACIPE core.topo -num_paras 7000 -DEID $i -DEFD 30 
		cd ..
	done
