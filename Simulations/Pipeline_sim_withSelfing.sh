  


###run sim before pipeline TBE


###transform all files 
#SC
bash ~/scripts/ParaMask_simulation_scripts/transform_mutation_SC.sh /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#SV
bash ~/scripts/ParaMask_simulation_scripts/transform_mutation_SV.sh /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2

####pipeline for single run
##convert to vcf
#SC
bash ~/scripts/ParaMask_simulation_scripts/run_genotypetype_to_mat_SC_withSelfing.sh 1 /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#SV
bash ~/scripts/ParaMask_simulation_scripts/run_genotypetype_to_mat_SV_new_withSelfing.sh 1 /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#concat SV duplicates
bash ~/scripts/ParaMask_simulation_scripts/concat_SV_mat_new.sh 1 /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#concat SV with SC
bash ~/scripts/ParaMask_simulation_scripts/concat_SC_SV_mat.sh 1 /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#annotate SVtab
bash ~/scripts/ParaMask_simulation_scripts/Annotate_SVtab.sh /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_final/SVtab_0.5_rep1.txt
#update postion with
bash ~/scripts/ParaMask_simulation_scripts/update_pos_sim_mat.sh SVtab_0.5_rep1_annotated.txt Simulations_SV_SC_rep1_mat.txt /netscratch/dep_coupland/grp_fulgione/bastiaan/data/SeDuS5_0.5_reps_withSelfing2
#simulate allelic ratios and coverage
bash ~/scripts/ParaMask_simulation_scripts/simulate_coverage_allele_depth.sh Simulations_SV_SC_rep1_mat_cpos.txt

