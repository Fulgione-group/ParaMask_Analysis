#!/bin/bash
awk 'BEGIN{FS=OFS="\t";c_sc=1; c_par=1}{if(FNR==1){$4="block"}else{if($3=="sc"){$4="SC_b"c_sc; c_sc++}else{$4="SV_b"c_par; c_par++}};print $0}' $1 > $(echo $1 | rev | cut -d'.' -f2- | rev)_annotated.txt
