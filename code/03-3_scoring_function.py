#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import sys
import argparse
import os
from scipy.stats import pearsonr, mannwhitneyu
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''
 
Description
    #########################################################
    This script calculate effect score of each DS-Rex_TX pair.
    The TPM file's index name MUST be "ID".
    OUTPUT: Whole_DS_score_{NMD,DOA,CDS_alt,and Whole}.txt
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat results', 
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument('--tu', '-t', 
                        help='Set whether using TU or not, default: Y',
                        required=False,
                        type=str,
                        default="Y")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

# %%
def Ref_cds(df):
    """ In case of Recovered, change their Sim TX CDS to Ref TX CDS
    """
    if df["pNMD"] != "-1":
        ref_cds = float(df["Stop"].split("-")[0]) - float(df["Start"].split("-")[0])
    else:   # if pNMD is Recovered, Ref CDS is calculated on the Sim TX
        ref_cds = float(df["Stop"].split("-")[1]) - float(df["Start"].split("-")[1])
    
    return ref_cds


def Direction(df):
    """ Assign directions for 2 possible classes

    Args:
        df (_type_): dataframe

    Returns:
        _type_: Inclusion case has 1 direction (same direction with dPSI) 
                other case has -1 direction (opposite direction with dPSI)
    """
    if df["event"] == "EI" or \
       df["event"] == "RI" or \
       df["event"] == "Alt_A3SS" or \
       df["event"] == "Alt_A5SS" or \
       df["event"] == "MXE2":
        return 1
    else:   # EI, RI, Alt3/5SS, MXE1
        return -1


def Make_int(input_df):
    """ Merge rMATS and interpretation to calculate mPSI

    Returns:
        _type_: merged_df
    """
    rmat = pd.read_csv("/".join(args.input.split("/")[:-2])+"/rmat.csv",
                       sep=",")
    rmat = rmat[["long_ID","FDR","dPSI_2_minus_1","mPSI1","mPSI2"]] # "mIC1","mSC1","mIC2","mSC2"
    min_fdr = rmat[rmat["FDR"] > 0]["FDR"].min()
    rmat["FDR"] = rmat["FDR"].apply(lambda x : x if x != 0 else min_fdr)
    rmat["-log10(q)"] = -np.log10(rmat["FDR"])
    rmat["PSI_direction"] = rmat["dPSI_2_minus_1"].apply(lambda x : 1 if x > 0 else -1)
    merged_df = pd.merge(input_df, rmat[["-log10(q)","long_ID","PSI_direction","dPSI_2_minus_1","mPSI1","mPSI2"]], left_on="LongID", right_on="long_ID")
    merged_df["domain_diff"] = merged_df["Ref_Domain"] - merged_df["Sim_Domain"]
    merged_df["Ref_TX_CDS"] = merged_df.apply(Ref_cds, axis=1)
    merged_df["ORF"] = merged_df["ORF"].apply(lambda x : 1 if x=="pORF1" else( 0.5 if x == "pORF2" else 0.1))
    merged_df["direction"] = merged_df.apply(Direction, axis=1)

    return merged_df


def Make_TPM():
    path = "/".join(args.input.split("/")[:-2])+"/tpm/" # TPM folder
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) if f.endswith(".tpm")]
    tu_check = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) if f.endswith("TU.txt")]

    if len(tu_check) > 0:
        return print("\n\nTU file is already existed\n\n")    

    elif len(files) > 0:  # In case of there are individual tpm file
        for k,file in enumerate(files):
            tpm = pd.read_csv(path+file,
                              sep="\t")
            tpm.index = tpm.iloc[:,0]
            tpm = tpm.iloc[:,1]
            if k == 0:
                merged_tpm = tpm
            else:
                merged_tpm = pd.merge(merged_tpm, tpm, 
                                      left_index=True, right_index=True,
                                      how="outer")
    
    else:   # In case of there is merged TPM file
        merged_tpm = pd.read_csv(path+"tpm_matrix.tsv",
                                 sep="\t")
        merged_tpm.index = merged_tpm.iloc[:,0]
        merged_tpm = merged_tpm.iloc[:,1:]
            
    return merged_tpm.fillna(0.0)


def Add_TU(tpm, query_df):
    """
    query_df = merged_input_for_eval.txt
    Calculate TU, using TPM value.
    First, calculate mean TPM for each transcript
    Then, using those values, TU is calculated in each gene
    """
    import time
    start_time = time.time()

    path = "/".join(args.input.split("/")[:-2])+"/tpm/" # TPM folder
    tu_check = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) if f.endswith("TU.txt")]
    if len(tu_check) > 0:   # TU file sholud be refrehsed
        return print("\n\nTU file is already existed\n\n")

    else:
        print("\n\nMaking TU...\n\n")
        tpm["mean"] = tpm.mean(axis=1)
        tpm = tpm[["mean"]]
        
        gene_tx = pd.read_csv("/".join(args.input.split("/")[:-2])+"tx_gene_dict",
                              sep="\t",
                              header=None,
                              index_col=0)
        if len(gene_tx[2].unique()) > 1:
            gene_tx = gene_tx[gene_tx[2].isin(query_df["Gene symbol"].unique())].copy()
        elif len(gene_tx[1].unique()) > 1:
            gene_tx = gene_tx[gene_tx[1].isin(query_df["Gene symbol"].unique())].copy()
            gene_tx.columns = [2,1] # In the subsequent steps, all script focus on col 2
        else:
            print("\nERROR: Check your tx_gene_dict file whether it has at least one gene annotation key\n")
            sys.exit(1)
            
        tu_input = pd.merge(tpm, gene_tx, 
                            left_index=True, 
                            right_index=True).sort_values(by=2).reset_index()
        print(f"Caldulated TU for {len(tu_input[2].unique())} genes")

        output_path = "/".join(args.input.split("/")[:-2]) + "/tpm/TU.txt"
        with open(output_path, 'w') as f:
            f.write("0"+"\t"+"1"+"\t"+"TU"+"\n")
            for gene in tu_input[2].unique():
                tx = tu_input[tu_input[2] == gene]
                tu_values = (tx["mean"] / tx["mean"].sum()).fillna(0.0).values
                output_data = np.column_stack((tx[["index",2]].values, tu_values))
                for row in output_data:
                    f.write("\t".join(map(str, row)) + "\n")
        
        end_time = time.time()
        print(f"Execution time (making TU): {(end_time - start_time) / 60} minutes")


def Multi_splicing_prop(score_df, stype, mode):
    """ Multiply splicing possibility

    Args:
        score_df (_type_): dataframe
        stype (_type_): NMD or DOA
        mode (_type_): w/wo TU (Y or N)
    """
    try:
        os.remove(args.input+"Classified_DS_{}.txt".format(stype))
        os.remove(args.input+"Whole_DS_score_{}.txt".format(stype))
    except:
        pass
    
    k = 0

    if mode == "N": # Make ones matrix
        tu = pd.DataFrame(np.ones(len(score_df["ref_tx"].unique())),
                          index=score_df["ref_tx"].unique(),
                          columns=["TU"])
    else:
        tu = pd.read_csv("/".join(args.input.split("/")[:-2])+"/tpm/TU.txt",
                         sep="\t",
                         index_col="0")

    for sid in score_df["LongID"].unique():    # Per each splicing
        sub_sid = score_df[score_df["LongID"]==sid] # Diverse simulated splicing consequence with matched isoform
        sub_sid = pd.merge(sub_sid, tu, left_on="ref_tx", right_index=True)
        if sub_sid.shape[0] > 0:
            for sub_class in sub_sid["event"].unique():    # Inclusion or Skipping
                each_SP = sub_sid[sub_sid["event"]==sub_class]

                if int(each_SP["direction"].unique()[0]) == 1:  # Inclusion cases
                    # Ps = ((each_SP["mPSI1"]+each_SP["mPSI2"])/2) # TPM * PSI
                    Ps = max([each_SP["mPSI1"].unique()[0], each_SP["mPSI2"].unique()[0]]) # Max PSI
                else:   # Skipping cases
                    # Ps = (((1 - each_SP["mPSI1"]) + (1 - each_SP["mPSI2"]))/2)
                    Ps = min([each_SP["mPSI1"].unique()[0], each_SP["mPSI2"].unique()[0]])  # Max (1-PSI)
                    Ps = 1 - Ps
    
                for tx in each_SP["ref_tx"].unique():   # Calculate Effect_Score per TX in each DS
                    ind_SP = each_SP[each_SP["ref_tx"]==tx].copy()
                    ind_SP["Effect_Score"] = (ind_SP["FS"].astype(float) * ind_SP["TU"] * Ps * abs(ind_SP["dPSI_2_minus_1"]) * ind_SP["ORF"].astype(float))
                    
                    if k == 0:
                        final_df = ind_SP
                    else:
                        final_df = pd.concat([final_df, ind_SP]) # Contains all case
                    k += 1

        else:   # No TPM isoform
            pass

    final_df = final_df.drop(["pair_ID", "ComplexID", "long_ID"], axis=1)
    final_df = final_df[["LongID","Gene symbol","ref_tx","event","Effect_Score","change_rate","pNMD","DOA_direction","dPSI_2_minus_1","TU","ORF","Start","Stop","dAA","5'UTR","3'UTR","Domain_integrity","Sim_Domain","Ref_Domain","-log10(q)"]]
    final_df.columns = ["LongID","Gene symbol","Reference_transcript","Simulated_event","Effect_Score","Domain_change_rate","Probability_of_NMD","Functional_class","Delta_PSI","Transcript_usage","ORF","AUG (Ref-Sim)","Stop_codon (Ref-Sim)","CDS_difference (nt)","5'UTR_difference (nt)","3'UTR_difference (nt)","Domain_integrity","Length_of_simulated_tx_domain","Length_of_referece_tx_domain","rMATS_FDR(-log10)"]
    final_df.to_csv(args.input+"Effect_score.tsv".format(stype),
                    sep="\t",
                    index=False)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/result/"):
        args.input = args.input+"/result/"
    
    ## Run code
    ref_result = pd.read_csv(args.input+"merged_input_for_eval.txt",
                             sep="\t")  # From 03-2_overview_change_rate.py
    merged_df = Make_int(ref_result)

# %%
    ## To calculate whole distribution of change rate, it contains NMD and DOA
    merged_df = merged_df.dropna(subset="Domain_integrity")
    merged_df = merged_df[merged_df["Domain_integrity"]!="Frame loss"]
    merged_df["change_rate"] = merged_df["change_rate"].astype(float)
    merged_df["change_rate"] = merged_df["change_rate"].apply(lambda x : 1 if x > 1 else x)
    
    ## Normalization of ∆CDS
    sub_data = merged_df[["Start","Stop","dAA"]]
    sub_data["CDS_ref"] = sub_data["Stop"].str.split("-").str[0].astype(float)-sub_data["Start"].str.split("-").str[0].astype(float)
    sub_data["CDS_sim"] = sub_data["Stop"].str.split("-").str[1].astype(float)-sub_data["Start"].str.split("-").str[1].astype(float)
    sub_data["dCDS"] = sub_data["CDS_ref"] - sub_data["CDS_sim"]
    
    def norm_CDS(df):
        if df["dCDS"] < 0:
            nCDS = abs(df["dCDS"]) / df["CDS_sim"]
        else:
            nCDS = abs(df["dCDS"]) / df["CDS_ref"]
    
    
        return nCDS
    sub_data["norm_CDS"] = sub_data.apply(norm_CDS,
                                          axis=1)
    merged_df["norm_CDS"] = sub_data["norm_CDS"]
    merged_df["FS"] = merged_df["change_rate"].astype(float) + merged_df["norm_CDS"].astype(float)

    ## Make TPM matrix
    if args.tu in ["Y","Y_own"]:
        merged_tpm = Make_TPM()
        Add_TU(merged_tpm, merged_df)

    Multi_splicing_prop(merged_df,"Whole",args.tu) 
