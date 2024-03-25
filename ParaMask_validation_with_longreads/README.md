# ParaMask validation with longreads

To validate ParaMask calls using long reads, our approach involved two key steps:

Firstly, we concentrated on assessing the overlap between ParaMask calls and duplications identified by long-read structural variation (SV) calling tools. This evaluation encompassed both multicopy and single-copy regions within the genome.

Secondly, we conducted a comparison between ParaMask calls and repeats identified from whole genome assemblies, incorporating not duplications only but all structural variants identified. This assessment spanned SV calling as well as assembly processes, utilizing PacBio HiFi long-read data from two distinct accessions belonging to the respective Arabis alpina population used to generate ParaMask calls.


This directory contains all scripts required to validate ParaMask calls using long-read data in accordance with our workflow.

## Assembly SV calling
To commence the validation workflow, we generated the whole genome assemblies and structural variant (SV) calls. The assembly process began with heterozygosity estimation using Jellyfish, followed by the creation of two independent assemblies using Hifiasm and Flye. Subsequently, the assemblies underwent purification via PurgeDups, followed by merging using Quickmerge. Finally, the resulting contigs were anchored against the reference chromosome to produce the final assembly. 
The SV calling process included mapping raw long reads against the reference genome, followed by the independent generation of two SV sets utilizing the tools cuteSV and sniffles. Afterwards, both SV sets underwent merging using SURVIVOR to yield a high-confidence SV set. Each step in these processes is represented by a corresponding script named after the program utilized and the specific step involved.

## Coverage Analysis longread assemblies




The directory Coverage Analysis longread assemblies contains all scripts used to annotate high coverage peaks inferred based on the alignment between assembled scaffolds and their corresponding longreads. Additionally, files of the RepeatMasker Analysis and the results of the Blast Search of the extracted high coverage peaks were used as input for the annotation. 

## Fit coordinate systems
This directory hosts a script designed for mapping scaffolds from generated assemblies to a reference. The script facilitates aligning coordinates between the assembly and the reference by utilizing the resultant mapping file. 

## Agreement analysis ParaMask SV
The directory Agreement analysis contains all scripts to validate the overlap between ParaMask and SV calls seperated for single-copy and multicopy regions. Each subdirectory contains a R-script to define a function how to calculate the agreement in the particular region. The single-copy directory additionally contains R-scripts to filter the ParaMask calls for ES03 and ES04 seed calls only. 

## RepeatMasker
To detect repeat elements in the generated assemblies, first the script RepeatMasker.sh was used. These repeats were afterwards lifted over to the reference using the script Bet_to_gff.R and Liftover_repeats.sh, which are all included in the directory RepeatMasker. 

## Agreement analysis ParaMask SV 

## Repeat SV overlap analysis

