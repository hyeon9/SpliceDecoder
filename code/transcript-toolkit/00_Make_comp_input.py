# %%
import pandas as pd
import numpy as np
import matplotlib
import subprocess
import argparse
import os
import time
import sys


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
start_time = time.time()

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''

Description
    #########################################################
    This script makes simulated transcript (Sim_TX) by using 
    given splicing events on the Ref_TX.
    #########################################################''',

    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)
    parser.add_argument("--cpat", "-cp", help="cpat path", 
                        required=True,
                        type=str)
    parser.add_argument("--cpatdb", "-cpdb", help="cpat db path", 
                        required=True,
                        type=str)
    parser.add_argument("--bedtools", "-b", help="bedtools path", 
                        required=True,
                        type=str)
    parser.add_argument("--genomefa", "-f", help="reference fasta", 
                        required=True,
                        type=str)
    parser.add_argument("--species", "-p", help="species name", 
                        required=True,
                        type=str)
    parser.add_argument('--query_list', '-q', 
                        required=True,
                        type=str)
        
    ## Optional
    parser.add_argument("--threads", "-t", help="number of threads to use, default: 9", 
                        type=int,
                        default="9")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


args, parser = parse_args(sys.argv[1:])
args.species = args.species.capitalize()
exon_dict = pd.read_csv(args.input+"exon_only.final.gtf",
                       sep="\t")    # From Final/00-2_processing_gtf.py
target_list = pd.read_csv(args.query_list,
                          sep="\t",
                          header=None)
exon_dict = exon_dict[exon_dict["ENSTID"].isin(target_list[0])]

cmd = "rm -r " + args.input + "cpat/"
subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
cmd = "rm -r " + args.input + "sim_bed/"
subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

# %%
def Make_output(paths):
    if not os.path.exists(paths):
        os.mkdir(paths)

sim_bed = f"{args.input}sim_bed/"
fasta_DIR = f"{args.input}get_Fasta/"
cpat_DIR = f"{args.input}cpat/"
Make_output(sim_bed)
Make_output(fasta_DIR)
Make_output(cpat_DIR)


def Run(tx_list, key):
    def Make_bed(tx):
        """ Make a bed file for the given transcripts (tx) in the given GTF file

        Args:
            tx (_type_): Transcript ID
        """
        
        df = exon_dict[exon_dict["ENSTID"]==tx]
        df = df[["chr","start","end","ENSGID","gene_symbol","strand","ENSTID","exon_number"]]
        df = df[df["ENSTID"]==tx].sort_values(by="exon_number")
        
        ## REMOVE
        # if df["strand"].unique()[0] == "-": #
        #     df = df.sort_values(by="start", ascending=False)
        #     df["exon_number"] = [i+1 for i in range(df.shape[0])]
        # else:
        #     pass
        
        df.loc[:,"start"] = df.loc[:,"start"].astype(int)-1
        df.to_csv(f"{sim_bed}{tx}.bed",
                  sep="\t", index=None,
                  header=None)
        # df.to_csv(f"{sim_bed}merged.bed",
        #           sep="\t", index=None,
        #           header=None,
        #           mode="a")


    def Get_seq(tx):
        """ Create Fasta file from the bed file (Make_bed)
        
        Args:
            tx : Transcript ID
        """
        
        if not os.path.exists(fasta_DIR):
            os.mkdir(fasta_DIR)

        ## Make AA sequence
        cmd = args.bedtools + " getfasta -fi " + args.genomefa + \
            " -bed " + sim_bed + tx + ".bed" + \
            " -fo " + fasta_DIR + tx + ".fa" + \
            " -s"
        subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

        files = open(fasta_DIR + tx + ".fa", "r")
        cpat_input = open(fasta_DIR + tx + "_cpat.fa", "w")
        merge_seq = ""
        for line in files:
            if not line.startswith(">"):
                merge_seq = merge_seq + line[:-1]
        
        proc_merged_seq = merge_seq.upper()
        cpat_input.write(">"+tx)
        cpat_input.write("\n"+proc_merged_seq)
        cpat_input.close()


    def ORF(tx):
        """ Predict ORF by using CPAT and the created Fasta (From Get_seq)
        
        Args:
            tx : Transcript ID
        """
        
        if not os.path.exists(cpat_DIR):
            os.mkdir(cpat_DIR)
        n_orf = 1
            
        ## Find ORF
        cmd = args.cpat + " -x " + args.cpatdb + args.species + "_Hexamer.tsv" + \
                " -d " + args.cpatdb + args.species + "*_logit*RData" + \
                " --verbose 0" + " --min-orf 1" + " --top-orf={}".format(n_orf) + \
                " -g "+ fasta_DIR + tx + "_cpat.fa" + \
                " -o" + cpat_DIR + tx + ".orf"
        subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

        try:    # Best ORF does not exist, fill zeros
            no_orf = pd.read_csv(f"{cpat_DIR}{tx}.orf.no_ORF.txt",
                                    sep="\t")  # If no_ORF exists, output will be filled by zeros
            cpat_output = open(f"{cpat_DIR}{tx}_bestorf.tsv", "a")
            cpat_output.write(tx+"\t"+
                            str(0)+"\t"+
                            str(0)+"\t"+
                            str(0)+"\t"+
                            str(0)+"\t"+
                            str(0)+"\n")   # Add matched tx, nt length of frame, coding prob, eventid
            
            cpat_output.close()


        except: # Best ORF exists
            files = open(f"{cpat_DIR}{tx}.orf.ORF_prob.tsv", "r")
            cpat_output = open(f"{cpat_DIR}{tx}_bestorf.tsv", "a")
            for line in files:
                if line.split("\t")[9].startswith("Coding_prob"):
                    pass
                else:
                    cpat_output.write(tx+"\t"+
                                      line.split("\t")[1]+"\t"+
                                      line.split("\t")[4]+"\t"+
                                      line.split("\t")[5]+"\t"+
                                      line.split("\t")[6]+"\t"+
                                      line.split("\t")[9])   # 1: tx length, 4,5: ORF start/end, 6: ORF length, 9: P coding
            cpat_output.close()            
    
    for tx in tx_list["ENSTID"].unique():    
        Make_bed(tx)
        Get_seq(tx)
        ORF(tx)


# %%
import time, os
from joblib import Parallel, delayed
n_cpu = args.threads
df_split = np.array_split(exon_dict, n_cpu) # Split into n dataframe
result = Parallel(n_jobs=n_cpu)(delayed(Run)(sub_df, str(k)) for k,sub_df in enumerate(df_split))
# %%
print("{} seconds".format(time.time() - start_time))
# %%
