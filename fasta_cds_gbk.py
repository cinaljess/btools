"""
This module looks at genbank features for coding sequences 'CDS' and then
outputs the gene name, orientation, and sequence in nucleotides per line. 
Infile searching is done utilizing the Bio python imports SeqIO and Seq.

First look up all features with Biopython, finding CDS, then get the location.
The orientation will determine the sequence as pulled from the below fasta and,
if complementary, the sequence will be complemented.

python fasta_cds_genbank.py GU19504.gb > genes_seq.txt
"""

import sys
import re
from Bio import SeqIO, Seq 

dat = SeqIO.read(open(sys.argv[1]), "genbank") 
genes = {}
for feat in dat.features:
    if feat.type == 'CDS':
        for line in str(feat).split("\n"):
            if "location" in line:
                loc = re.findall("[0-9]*:[0-9]*", line)
                loc = [x for x in loc if len(x) > 1][0]
                ori = re.findall("\(.\)", line)[0].strip("(").rstrip(")")
        pseq = feat.qualifiers['translation'][0]
        genes[feat.qualifiers['gene'][0]] = [pseq, loc, ori]

for gene in genes:
    seq = ''
    if genes[gene][2] == "+":
        start, end = genes[gene][1].split(":")
        seq = dat.seq[int(start):int(end)]
    else:
        start, end = genes[gene][1].split(":")
        seq = dat.seq[int(start):int(end)]
        seq = Seq.reverse_complement(seq)
    genes[gene].append(seq)
    print(gene, genes[gene][2], str(genes[gene][3]))
