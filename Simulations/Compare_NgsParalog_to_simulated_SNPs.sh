##### SCRIPT TO GENERATE SIMULATED ILLUMINA PAIRED END READS ###### 30.07.2025
## ====== 1) Create Inputfile with all sample name in col1 and sample name + hap ID in col2: ======== ##
wrdir="/path/to/Paramask_sim_fastas/"
cd ${wrdir}
for f in *_hap[12].fasta; do
  hap=${f%.fasta}
  base=${hap%_hap[12]}
  echo -e "${base}\t${hap}"
done > ${wrdir}ALL_samples.txt


##### =========== FIRST GENERATE FOR EACH SAMPLE A FILE WITH ALL IDs INSIDE ============== ########
coverage=5      # Adjust for coverage 10X
samples=15      # Adjust for sample size 100
condition_vector=("HW_Rep1" "HW_Rep2" "HW_Rep3" "Fis0.9_Rep1" "Fis0.9_Rep2" "Fis0.9_Rep3")
for condition in ${condition_vector[@]}; do
  input_IDs="${wrdir}${condition}_${samples}N_samples.txt"
  grep "$condition" ${wrdir}ALL_samples.txt | shuf | head -n $samples > $input_IDs
done



### ======== FROM HERE FINAL SCRIPT UNTIL SAMTOOLS PILEUP ============= ###
coverage=5
samples=15
condition_vector=("HW_Rep1" "HW_Rep2" "HW_Rep3" "Fis0.9_Rep1" "Fis0.9_Rep2" "Fis0.9_Rep3")
for condition in ${condition_vector[@]}; do
  ## ======= 1) generate fastq for each of the samples ============== ##
  echo "generate fastq for each of the samples"
  input_IDs="/path/to/Paramask_sim_fastas/${condition}_${samples}N_samples.txt"
  ## WAIT UNTIL 20 JOBS are DONE:
  max_jobs=20
  job_count=0

  while read -r col1 col2; do
    wrdir="/path/to/Paramask_sim_fastas"
    base_path_outdir="/path/to/ParaMask/Seq_sim/${coverage}X_${samples}N/${condition}"
    out_dir="${base_path_outdir}/${col1}"
    mkdir -p "${out_dir}"

    (
      cd "${out_dir}" || exit
      outfile_prefix="${col2}"
      sequencing_system="HSXt"
      seq_ref_file="${wrdir}/${col2}.fasta"

      echo "Simulating ${col2} in ${out_dir}"
      LD_LIBRARY_PATH=$HOME/lib art_illumina -ss "$sequencing_system" -p -i "$seq_ref_file" \
        -l 150 -f ${coverage} --sdev 10 --noALN -o "$outfile_prefix" --mflen 200
    ) &

    job_count=$((job_count + 1))
    if [[ $job_count -ge $max_jobs ]]; then
      wait   # wait bis 20 jobs fertig sind
      job_count=0
    fi
  done < $input_IDs

  wait


  ## ======= 2) merge all the fastq files together into one  ============== ##
  echo "merge all the fastq files together into one"
  input_IDs="/path/to/Paramask_sim_fastas/${condition}_${samples}N_samples.txt"
  base_path_outdir="/path/to/ParaMask/Seq_sim/${coverage}X_${samples}N/${condition}"
  ngsparalog_output="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/bams"
  reference_fasta="/path/to/Paramask_sim_fastas/reference.fasta"

  mkdir -p "$ngsparalog_output"
  mkdir -p "$base_path_outdir"

  job_count=0
  max_jobs=20

  awk '{print $1}' "$input_IDs" | sort | uniq | while read -r col1; do

    (
      out_dir="${base_path_outdir}/${col1}"
      echo "Merging FASTQs for ${col1} in ${out_dir}"
      cat "${out_dir}"/*1.fq > "${out_dir}/${col1}_merged_1.fq"
      cat "${out_dir}"/*2.fq > "${out_dir}/${col1}_merged_2.fq"

      forward="${out_dir}/${col1}_merged_1.fq"
      reverse="${out_dir}/${col1}_merged_2.fq"
      echo "Mapping of ${col1}"
      bwa mem -M -t 8 "$reference_fasta" "$forward" "$reverse" | samtools view -b - | samtools sort -@ 8 -o "${ngsparalog_output}/${col1}.bam"
    ) &

    job_count=$((job_count + 1))
    if [[ $job_count -ge $max_jobs ]]; then
      wait
      job_count=0
    fi

  done
  wait  # Ensure all merges and mappings are complete

  ## ========== 3) generate input for ngsparalog: bam list und position file ================ ##
  echo "generate input for ngsparalog: bam list und position file"

  bam_list="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}_bam_list.txt"
  path="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/bams"

  cd "$path"
  for f in *.bam; do
    bam="${path}/${f}"
    echo "${bam}"
  done > "$bam_list"

  file="/path/to/Paramask_sim_fastas/Position_files_${condition}.txt"
  position_list="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}_positions.txt"
  awk '{print $1}' "$file" | while read -r col1; do
    echo -e "reference\t${col1}"
  done > "$position_list"
  wait

  ## ======= 4) input the position file into samtools mpileup to generate pileup file ========== ##
  echo "input the position file into samtools mpileup to generate pileup file"
  ngs_in_file="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}.pileup"

  (
    samtools mpileup -b "$bam_list" -l "$position_list" --ff UNMAP,DUP > "$ngs_in_file"
  ) &

  job_count=$((job_count + 1))
  if [[ $job_count -ge $max_jobs ]]; then
    wait
    job_count=0
  fi

  wait

  ## ========== 5) Use piöeup file to run NgsParalog: =========== ##
  for condition in ${condition_vector[@]}; do
    position_list="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}_positions.txt"
    bam_list="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}_bam_list.txt"
    ngs_in_file="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}.pileup"
    ngs_out_file="/path/to/ParaMask/ngsParalog/seq_sim/${coverage}X_${samples}N/${condition}/${condition}.lr"
    /path/to/software/ngsParalog/ngsParalog calcLR -infile ${ngs_in_file} -outfile ${ngs_out_file} -minQ 20 -minind 5 -mincov 1 &
  done
  wait

  echo "All files of NgsParalog are generated"

  ## Continue in R with processing using the Script: Processing_NgsParalog_sim_comparison.R ##
done
