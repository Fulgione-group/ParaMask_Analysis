#!/bin/bash

query="/Path/to/query_assembly.fasta"
reference="/Path/to/reference_assembly.fasta"
quickmerge="/Path/to/quickmerge/"


export PATH=${quickmerge}:$PATH
merge_wrapper.py ${reference} ${query}
