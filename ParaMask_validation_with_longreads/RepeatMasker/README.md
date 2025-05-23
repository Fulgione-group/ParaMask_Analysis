# Overview

This directory contains all scripts necessary for detecting repeat elements in an assembly and subsequently transferring those repeats onto another assembly to align coordinates between the two assemblies.

The script _RepeatMasker.sh_ [Smit et al., 2013][https://github.com/rmhubley/RepeatMasker] is employed to annotate repeats based on an input assembly.

Following this, the script _Bed_to_gff.R_ transforms the resulting repeat file into a GFF file, which serves as input for the script _Liftover_repeats.sh_ [Shumate et al., 2021][https://github.com/agshumate/Liftoff]. This latter script facilitates the transfer of repeats from one assembly to another, thereby inferring repeat coordinates of the target assembly."
