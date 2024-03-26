  
###transform all files 
#SC
bash transform_mutation_SC.sh PATH_TO_OUT_DIR
#SV
bash transform_mutation_SV.sh PATH_TO_OUT_DIR

####pipeline for single run
##convert to vcf
#SC
bash run_genotypetype_to_mat_SC_withSelfing.sh 1 PATH_TO_OUT_DIR
#SV
bash run_genotypetype_to_mat_SV_new_withSelfing.sh 1 PATH_TO_OUT_DIR
#concat SV duplicates
bash concat_SV_mat_new.sh REP_NUMBER PATH_TO_OUT_DIR
#concat SV with SC
bash concat_SC_SV_mat.sh REP_NUMBER PATH_TO_OUT_DIR
#annotate SVtab
bash Annotate_SVtab.sh SVtab.txt
#update postion with
bash update_pos_sim_mat.sh SVtab_annotated.txt Simulations_SV_SC_rep{REP_NUMBER}_mat.txt PATH_TO_OUT_DIR 
#simulate allelic ratios and coverage
bash ~/scripts/ParaMask_simulation_scripts/simulate_coverage_allele_depth.sh Simulations_SV_SC_rep{REP_NUMBER}_mat_cpos.txt

