#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
from fasta_checker import Fasta_checker

db = "16sRNA"
query = "../../../Projects/data/genomes/2018-11-06-Coxiella_burnetii_RSA_493/input/GCF_000007765.2_ASM776v2_genomic.fa"
query2 = "../../../Projects/data/genomes/2019-06-06-Coxiella_burnetii_California_16_RSA350_clone_2/input/GCF_002924325.1_ASM292432v1_genomic.fa"
fc = Fasta_checker(query2, db)
#blast_out = fc.blast()
if fc.check_contigency():
    print("cont")
else:
    print("not_cont")
fc.parse_blast(db)
