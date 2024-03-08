#!/bin/bash

###
#		Get JSFS
###

# alpina paramask:
matrix=insert_path_to_SNP_matrix
resultsDir=insert_path_to_results_folder
mkdir -p ${resultsDir}

for samples1 in insert_path_to_samples_file1
do
	for samples2 in /insert_path_to_samples_file2
	do
		samplesOutgroup=insert_path_to_outgroup_file
		for mask in /insert_path_to_mask_file
		do
			results=${resultsDir}sfs/
			mkdir -p ${results}
			echo ${results}
			
			java -Xmx4G -classpath ~/java/lib/junit.jar:~/java/lib/jbzip2-0.9.1.jar:~/java/lib/args4j-2.0.12.jar:~/java/lib/commons-compress-1.0.jar:~/java/lib/gson-1.6-javadoc.jar:~/java/lib/gson-1.6-sources.jar:~/java/lib/gson-1.6.jar:bin c.e.data_processing.Matrix_getJSFS_general_2popsPlusOutgroup ${matrix} ${mask} ${results}${short1}-${short2} ${samples1} ${samples2} ${samplesOutgroup}
		done
	done
done
