#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-
from Bio.Blast.Applications import NcbiblastnCommandline
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO
import sys


class Fasta_checker():
    def __init__(self, fasta, gene):
        out_path = "./out"
        db_path = "./db"
        self.query = fasta
        self.db = os.path.join(db_path, "{}.fasta".format(gene))
        self.blast_out = os.path.join(out_path, "{}.xml".format(gene))

    def check_sequence_number(self):
        sn = len([len(rec) for rec in SeqIO.parse(self.query, "fasta")])
        if sn > 2:
            sys.stdout.write("{} contains more than 2 chromosomes, it'S Most likely a contig \n".format(self.query))
            sys.stdout.flush()
        else:
            sys.stdout.write("{} contains 1 chromosome, it'S Most likely a complete genome \n".format(self.query))
            sys.stdout.flush()
            

    def blast(self):
        blastx_cline = NcbiblastnCommandline(query=self.query, db=self.db,
                                             evalue=0.001,
                                             outfmt=5, out=self.blast_out)
        blastx_cline()
        return self.blast_out

    def parse_blast(self, db):
        print("..........checking {} sequence...........".format(db))
        if os.path.exists(self.blast_out):
            print(self.blast_out)
            with open(self.blast_out) as xml_out:
                blast_record = NCBIXML.parse(xml_out)
                for record in blast_record:
                    for items in record.alignments:
                        for hsp in items.hsps:
                            if hsp.expect < 0.0:
                                if db == "16sRNA":
                                    expected_coverage = 1456
                                    if hsp.identities == expected_coverage:
                                        print("{} check was succesful")
                                    else:
                                        print("expected coverage is {} \n but the observed coverage is {}".format(expected_coverage, hsp.identities)
                                              )
                  #  for hsp in alignment.hsps:
                  #      print(hsp.expect)

    def check_contigency(self):
        contig_order = True
        contigency_list = []
        nsc = 0
        for record in SeqIO.parse(self.query, "fasta"):
            if nsc < 10:
                contigency_list.append(len(record.seq))
                nsc += 1
            else:
                break
        print(contigency_list)
        assumed_max = contigency_list[0]
        for size in contigency_list:
            if size > assumed_max:
                return False
        return contig_order
