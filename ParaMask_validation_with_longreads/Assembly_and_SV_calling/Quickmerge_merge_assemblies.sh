#!/bin/bash

query="query_assembly"
reference="reference_assembly"

export PATH=/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge:/netscratch/dep_coupland/grp_fulgione/male/assemblies/it20_007/pacbio_hifi/merge_hifiasm_flye/quickmerge/MUMmer3.23:$PATH
merge_wrapper.py ${reference} ${query}
