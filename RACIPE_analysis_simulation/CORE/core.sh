#!/bin/bash
#cd CORE
#./RACIPE *.topo -num_paras 7000 &
for i in {2..3}
	do 
		mkdir set_$i
		cd set_$i
		scp ../../RACIPE_dummy/* .
		scp ../../core.topo .
		./RACIPE core.topo -num_paras 7000 &> core_$i.out & 
		cd ..
	done
