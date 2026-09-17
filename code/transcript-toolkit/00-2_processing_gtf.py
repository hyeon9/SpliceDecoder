#!/usr/bin/env python3
# %%
import argparse
import sys

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
    description='''
Description
    #########################################################
    Raw gencode gtf should be processed using below command before run this script
    "awk -F "\t" '{if ($3 == "exon") print }' ${Your.gtf} > exon_only.gtf"
    This script makes proper format of exon gtf
    #########################################################
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--exon', '-e', 
                        help='Processed exon gtf file, it should contains all exon number', 
                        required=True,
                        type=str)
    parser.add_argument('--out', '-o', 
                        help='Output directory', 
                        required=True,
                        type=str)

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

def Make_gtf(input,dir):
    """ Processing exon only gtf file
        Raw gencode gtf should be processed in shell script before run this script
        "awk -F "\t" '{if ($3 == "exon") print }' gencode.v32.basic.annotation.gtf > exon_only.gtf"
    """

    gtf = open(input)
    out = open(dir+"/exon_only.final.gtf", "w")
    out.write("chr""\t""start""\t""end""\t""strand""\t""ENSGID""\t""gene_symbol""\t""ENSTID""\t""exon_number""\n")
    for k,line in enumerate(gtf):
        if line.startswith("#"):
            pass
        else:
            line = line.strip()[:-1]
            header_check = []
            chr = str(line.split("\t")[0])
            st = str(line.split("\t")[3])
            ed = str(line.split("\t")[4])
            strand = str(line.split("\t")[6])
            if len(line.split("\t")[8].split(";")) < 4:   # Long read gtf (PacBio)
                for each in line.split("\t")[8].split(";"):
                    if each.strip().startswith("gene_id"):
                        g_id = str(each.strip().split("\"")[1])
                        g_symbol = str(each.strip().split("\"")[1])
                        header_check.append(g_id)
                    elif each.strip().startswith("transcript_id"):
                        tx_id = str(each.strip().split("\"")[1])
                        header_check.append(tx_id)
                    elif each.strip().startswith("exon_number"):
                        en = str(each.strip().split(" ")[1])
                        header_check.append(en)
                    else:
                        pass
            else:   # General gencode gtf
                for each in line.split("\t")[8].split(";"):
                    if each.strip().startswith("gene_id"):
                        g_id = str(each.strip().split("\"")[1])
                        g_symbol = str(each.strip().split("\"")[1])
                        header_check.append(g_id)
                    elif each.strip().startswith("transcript_id"):
                        tx_id = str(each.strip().split("\"")[1])
                        header_check.append(tx_id)
                    elif each.strip().startswith("exon_number"):
                        en = str(each.strip().split(" ")[1])
                        header_check.append(en)
                    else:
                        pass

            out.write(chr+'\t'+st+'\t'+ed+'\t'+strand+'\t'+g_id+'\t'+g_symbol+'\t'+tx_id+'\t'+en+
                      '\n')
    
    out.close()

def Make_gene_dict(DIR):
    import os
    files = [f for f in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, f)) if f.startswith("main") if f.endswith(".gtf") ]
    if len(files) == 0:
        print("Check your main.gtf")
    else:
        f = open(DIR+"{}".format(files[0]))
        o = open(DIR+"tx_gene_dict", "w")
        for line in f:
            read_line = line.strip().split("\t")
            if not line.startswith("#"):
                if read_line[2] == "transcript":
                    if "transcript_id" in read_line[8].split(" "):
                        tx = read_line[8].split(" ")[read_line[8].split(" ").index("transcript_id")+1].split("\"")[1]
                    if "gene_name" in read_line[8].split(" "):
                        gene = read_line[8].split(" ")[read_line[8].split(" ").index("gene_name")+1].split("\"")[1]
                    elif "gene_id" in read_line[8].split(" "):
                        gene = read_line[8].split(" ")[read_line[8].split(" ").index("gene_id")+1].split("\"")[1]
                    else:
                        gene = ""
                    o.write(tx+"\t"+gene+"\n")
        f.close()
        o.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    else:
        args, parser = parse_args(sys.argv[1:])    
        Make_gtf(args.exon, args.out)
        Make_gene_dict(args.out)

    gtf = open(input)
    out = open(dir+"/exon_only.final.gtf", "w")
    out.write("chr""\t""start""\t""end""\t""strand""\t""ENSGID""\t""gene_symbol""\t""ENSTID""\t""exon_number""\n")
    for k,line in enumerate(gtf):
        if line.startswith("#"):
            pass
        else:
            line = line.strip()[:-1]
            header_check = []
            chr = str(line.split("\t")[0])
            st = str(line.split("\t")[3])
            ed = str(line.split("\t")[4])
            strand = str(line.split("\t")[6])
            if len(line.split("\t")[8].split(";")) < 4:   # Long read gtf (PacBio)
                for each in line.split("\t")[8].split(";"):
                    if each.strip().startswith("gene_id"):
                        g_id = str(each.strip().split("\"")[1])
                        g_symbol = str(each.strip().split("\"")[1])
                        header_check.append(g_id)
                    elif each.strip().startswith("transcript_id"):
                        tx_id = str(each.strip().split("\"")[1])
                        header_check.append(tx_id)
                    elif each.strip().startswith("exon_number"):
                        en = str(each.strip().split(" ")[1])
                        header_check.append(en)
                    else:
                        pass
            else:   # General gencode gtf
                for each in line.split("\t")[8].split(";"):
                    if each.strip().startswith("gene_id"):
                        g_id = str(each.strip().split("\"")[1])
                        g_symbol = str(each.strip().split("\"")[1])
                        header_check.append(g_id)
                    elif each.strip().startswith("transcript_id"):
                        tx_id = str(each.strip().split("\"")[1])
                        header_check.append(tx_id)
                    elif each.strip().startswith("exon_number"):
                        en = str(each.strip().split(" ")[1])
                        header_check.append(en)
                    else:
                        pass

            out.write(chr+'\t'+st+'\t'+ed+'\t'+strand+'\t'+g_id+'\t'+g_symbol+'\t'+tx_id+'\t'+en+
                      '\n')
    
    out.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    else:
        args, parser = parse_args(sys.argv[1:])    
        Make_gtf(args.exon, args.out)
