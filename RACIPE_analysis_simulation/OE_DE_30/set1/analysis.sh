#!/bin/bash
for i in {1..6}
	do
	       cd DEid$i
       	       python3 ../phenotype_analysis.py
	       scp analysis/*hist_phenotype.png ../../hist_analysis_OE_DE_20_30_50
	       cd ..
        done	       
