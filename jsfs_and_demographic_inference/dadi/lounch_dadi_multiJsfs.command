#!/bin/bash
#

# python simple_unfolded.py 6 exp ./results_cvi_2020-04-21/

for sfs in path_to_jsfs
do
	results=${sfs}_results/
	mkdir -p ${results}
	
	for r in {0..199}
	do
  		python simple_unfolded.py ${r} 'pop_split' ${results} ${sfs} 
   		python simple_unfolded.py ${r} 'exp' ${results} ${sfs}
 		python simple_unfolded.py ${r} 'bottleneck' ${results} ${sfs}
 		python simple_unfolded.py ${r} 'im' ${results} ${sfs}
		python simple_unfolded.py ${r} 'bottleneck_twoSides' ${results} ${sfs}
	done
done

