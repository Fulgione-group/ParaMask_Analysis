This directory contains all scripts to detect repeat elemts in an assembly and to subsequently, lift those repeats over onto another assembly to fit coordinates between two assemblies. 

With the script _RepeatMasker.sh_, repeats are annotated based on an assembly as input. 

Using the script _Bed_to_gff.R_ the resulting repeat file is transformed to a gff file which can now be used for the script _Liftover_repeats.sh_. This script liftsover the repeats from one assembly on to the other, inferring repeat coordinates of the other assembly. 


This directory contains all scripts necessary for detecting repeat elements in an assembly and subsequently transferring those repeats onto another assembly to align coordinates between the two assemblies.

The script _RepeatMasker.sh_ is employed to annotate repeats based on an input assembly.

Following this, the script _Bed_to_gff.R_ transforms the resulting repeat file into a GFF file, which serves as input for the script _Liftover_repeats.sh_. This latter script facilitates the transfer of repeats from one assembly to another, thereby inferring repeat coordinates in the target assembly."
