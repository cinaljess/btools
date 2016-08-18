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
