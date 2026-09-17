#!/usr/bin/env python3
# %%
import argparse
import sys

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
        description='''
        
Description
    #########################################################
    This script adds exon number tag in your gtf.
    If your gtf file does not have exon number, Splice-decoder does not work properly
    #########################################################
    ''',
        formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains exon gtf file', 
                        required=True,
                        type=str)
    parser.add_argument('--seq', '-s', 
                        help='Sequencing type used in GTF, e.g., SR (short-read), LR (Long-read)', 
                        required=True,
                        type=str)
    
    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

## Assembled gtf (both LR and SR)
## awk -F "\t" '{if ($2 != "PacBio") print }' exon_only.gtf > non_pacbio_exon_only.gtf
## awk -F "\t" '{if ($2 == "PacBio") print }' exon_only.gtf > LR_exon_only.gtf

def Add_en(DIR,seq_type):
    if seq_type == "LR":
        f = open(DIR+"LR_exon_only.gtf")
    else:
        f = open(DIR+"exon_only.gtf")

    df = {}
    pre_tx = "start"
    en = 0
    df[pre_tx] = []
    
    ## Add exon number to tx in dict, using this dict in the below for loop
    for line in f:
        if pre_tx != line.split("\t")[-1][:-1]:
            df[pre_tx].append(en)
            df[line.split("\t")[-1][:-1]] = []
            en = 1
            pre_tx = line.split("\t")[-1][:-1]
        
        else:
            en += 1
    df[pre_tx].append(en)   # For end of the tx

    if seq_type == "LR":
        f = open(DIR+"LR_exon_only.gtf")
        o = open(DIR+"OUT_LR_exon_only.gtf", "w")

    else:
        f = open(DIR+"exon_only.gtf")
        o = open(DIR+"OUT_exon_only.gtf", "w")

    pre_tx = "start"
    for line in f:
        if pre_tx != line.split("\t")[-1][:-1]:
            if line.split("\t")[6] == "-":
                en = df[line.split("\t")[-1][:-1]][0]
                o.write("\t".join(line.strip().split("\t"))+" exon_number {};".format(str(en))+"\n")
            else:
                en = 1
                o.write("\t".join(line.strip().split("\t"))+" exon_number {};".format(str(en))+"\n")
            pre_tx = line.split("\t")[-1][:-1]
        else:
            if line.split("\t")[6] == "-":
                en -= 1
                o.write("\t".join(line.strip().split("\t"))+" exon_number {};".format(str(en))+"\n")
            else:
                en += 1
                o.write("\t".join(line.strip().split("\t"))+" exon_number {};".format(str(en))+"\n")
    o.close()

    import subprocess
    if seq_type == "LR":
        cmd = 'cat {}OUT_LR_exon_only.gtf \
               {}non_pacbio_exon_only.gtf \
               > {}exon_only.gtf'.format(DIR,DIR,DIR)
        subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    else:
        cmd = 'mv {}OUT_exon_only.gtf \
               {}exon_only.gtf'.format(DIR,DIR)
        subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

# %%
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"
    Add_en(args.input, args.seq)
