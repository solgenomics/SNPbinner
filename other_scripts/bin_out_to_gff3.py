#!/usr/bin/env python
"""
usage: bin_out_to_gff3.py [-h] [-s SEQ_ID] bin_file out

positional arguments:
  bin_file              Path to snnpbinnner bins output
  out                   Path for output file.

optional arguments:
  -h, --help            show this help message and exit
  -s SEQ_ID, --seq_id SEQ_ID
                        sequence ID for the gff file. Defaults to the binmap
                        ID.
"""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("bin_file", type=Path, help="Path to snnpbinnner bins output")
parser.add_argument("out", type=Path, help="Path for output file.")
parser.add_argument("-s","--seq_id", type=str, help="sequence ID for the gff file. Defaults to the binmap ID.", required=False, default="")
args = parser.parse_args()

seqid = args.seq_id
source = "snpbinner"
score = "."
strand = "."
phase = "."

attr_template = "Name={ril}.{num:04d}[{center}]"
gff_line = "{seqid}\t{source}\t{type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"

with args.bin_file.open() as input_csv, args.out.open("w") as output_gff:
    output_gff.write("##gff-version 3\n")
    line = "-"
    binmap_ids = []
    while not line.startswith("##bin start"): 
        if line.startswith("##binmap id"):
            binmap_ids = [index.strip() for index in line.split(",")[1:]]
        line = input_csv.readline()
    bin_starts = [int(index.strip()) for index in line.split(",")[1:]]
    while not line.startswith("##bin end"): line = input_csv.readline()
    bin_ends = [int(index.strip()) for index in line.split(",")[1:]]
    while not line.startswith("bin center"): line = input_csv.readline()
    bin_centers = [int(index.strip()) for index in line.split(",")[1:]]
    
    for line in input_csv:
        cells = line.split(",")
        if len(cells)<len(bin_starts)+1: continue
        ril = cells[0]
        binGenotypes = [geno.strip() for geno in cells[1:]]
        print(binGenotypes)
        for binNum,genotype in enumerate(binGenotypes):
            type = "genobin_"+genotype
            start = bin_starts[binNum]+1
            end = bin_ends[binNum]
            attributes = attr_template.format(ril=ril,
                                              num=binNum,
                                              center=bin_centers[binNum])
            out_line = gff_line.format(seqid=seqid if seqid!="" else binmap_ids[binNum],
                                       source=source,
                                       type=type,
                                       start=start,
                                       end=end,
                                       score=score,
                                       strand=strand,
                                       phase=phase,
                                       attributes=attributes)
            output_gff.write(out_line)
