#!/bin/bash
reference="/Path/to/reference.fasta"
assembly="/Path/to/assembly.fasta"
ID="Accession_ID"
Ragtag="/Path/to/ragtag/"
mapper="/Path/to/mapper/minimap"

export PATH="${Ragtag}/bin/:$PATH"

# without chimeric contig correction
ragtag.py scaffold ${reference} ${assembly} -o ${ID}_RagTag -C --aligner ${mapper} -t 20

