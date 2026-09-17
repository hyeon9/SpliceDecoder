# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import argparse
import os
import sys


def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''

Description
    #########################################################
    This script estimate functional difference between Ref_TX and Sim_TX
    This work can be focused on certain functional categoriy 
    (e.g., binding, dna_bind, motif, domain, region, and all)
    OUTPUT: _Domain_integrity_indi, _NMD_check, _Main_output
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)
    parser.add_argument('--query_list', '-q', 
                        required=True,
                        type=str)
    parser.add_argument('--config_path', '-p',
                        help='The path of your config file',
                        required=True,
                        type=str)
    parser.add_argument("--code_dir", "-c", help="Path of code",
                        required=True,
                        type=str)
    parser.add_argument("--species", "-s", help="Species",
                        required=True,
                        type=str)
    
    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args


def Load_data():
    sim_bed = args.input+"sim_bed/"
    cpat_DIR = args.input+"cpat/"

    ## coreID = SID|ENST|event_type
    if not os.path.exists(args.input+"table/"):
        os.mkdir(args.input+"table/")

    cmd = "cat "+cpat_DIR+"*bestorf.tsv > " + \
           args.input+"merged_bestorf.txt"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    ## Make domain data
    cmd = "cat "+sim_bed+"*bed > " + args.input + "merged.bed"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    cmd = "bash " + args.code_dir + "/run_intersect.sh " + args.config_path 
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    w_pfam = pd.read_csv(args.input + "table/Pfam.txt",
                       sep="\t", header=None)

    return w_pfam


def Make_query(final_bed, pfam):
    """ Collect genomic coordinate information from query event

    Args:
        final_bed (dataframe): 01-4 output (merged_bed)

    Returns:
        new_line : Input of the drawing figure function (in the relative position, start from zero)
    """
    diff = []
    new_line = []
    domain_line = []
    domain_name = []
    strand = final_bed[5].values[0]
    pfam = pfam[pfam[5]==strand]
    pfam[3] = pfam[3]+";"+pfam[4]+";"+pfam[6]
    for line in range(final_bed.shape[0]):
        if line > 0:
            if strand == "+":
                ## ADD strand specific domain assign and multiple domain case
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                intron = np.abs(int(final_bed.iloc[line-1,2] - int(final_bed.iloc[line,1])))
                diff.append(intron)
                new_line.append([int(final_bed.iloc[line,1])-np.sum(diff)-zero, 
                                 int(final_bed.iloc[line,2])-np.sum(diff)-zero])

                for n_domain in range(domain.shape[0]):
                    domain_line.append([int(domain.iloc[n_domain,1])-np.sum(diff)-zero,
                                        int(domain.iloc[n_domain,2])-np.sum(diff)-zero,])
                        
                    if pre_domain == domain.iloc[n_domain,3]:
                        domain_name.append(" ")
                        pre_domain = domain.iloc[n_domain,3]
                    else:
                        domain_name.append(domain.iloc[n_domain,3])
                        pre_domain = domain.iloc[n_domain,3]
        
            else:   # Start with -2, because first exon is -1
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                intron = np.abs(np.abs(int(final_bed.iloc[line-1,1])-zero) - np.abs(int(final_bed.iloc[line,2])-zero))
                diff.append(intron)
                new_line.append([np.abs(int(final_bed.iloc[line,2])-zero)-np.sum(diff), 
                                np.abs(int(final_bed.iloc[line,1])-zero)-np.sum(diff)])    # 2*zero, represents to zero and then direction change
                for n_domain in range(domain.shape[0]):
                    domain_line.append([np.abs(int(domain.iloc[n_domain,2])-zero)-np.sum(diff),
                                        np.abs(int(domain.iloc[n_domain,1])-zero)-np.sum(diff),])

                    if pre_domain == domain.iloc[n_domain,3]:
                        domain_name.append(" ")
                        pre_domain = domain.iloc[n_domain,3]
                    else:
                        domain_name.append(domain.iloc[n_domain,3])
                        pre_domain = domain.iloc[n_domain,3]

        else:   # First line
            pre_domain = " "
            if strand == "+":
                zero = int(final_bed.iloc[line,1])  # The most left exon's start site
                new_line.append([int(final_bed.iloc[line,1])-zero,
                                 int(final_bed.iloc[line,2])-zero])
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                
                for n_domain in range(domain.shape[0]):
                    domain_line.append([int(domain.iloc[n_domain,1])-zero,
                                        int(domain.iloc[n_domain,2])-zero,])
                    if n_domain == 0:
                        domain_name.append(domain.iloc[n_domain,3])
                        pre_domain = domain.iloc[n_domain,3]
                    else:
                        if pre_domain == domain.iloc[n_domain,3]:
                            domain_name.append(" ")
                            pre_domain = domain.iloc[n_domain,3]
                        else:
                            domain_name.append(domain.iloc[n_domain,3])
                            pre_domain = domain.iloc[n_domain,3]

            else:
                zero = int(final_bed.iloc[line,2])  # The most right exon's end site
                new_line.append([np.abs(int(final_bed.iloc[line,2])-zero),
                                 np.abs(int(final_bed.iloc[line,1])-zero)])
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3]].drop_duplicates()
                
                for n_domain in range(domain.shape[0]):
                    domain_line.append([np.abs(int(domain.iloc[n_domain,2])-zero),
                                        np.abs(int(domain.iloc[n_domain,1])-zero),])
                    if n_domain == 0:
                        domain_name.append(domain.iloc[n_domain,3])
                        pre_domain = domain.iloc[n_domain,3]
                    else:
                        if pre_domain == domain.iloc[n_domain,3]:
                            domain_name.append(" ")
                            pre_domain = domain.iloc[n_domain,3]
                        else:
                            domain_name.append(domain.iloc[n_domain,3])
                            pre_domain = domain.iloc[n_domain,3]
                
    return new_line, domain_line, domain_name


####### Debug
#from types import SimpleNamespace

#args = SimpleNamespace(
#    input="/home/kangh/lab-server/SpliceDecoder/Isoformswitch/",
#    config_path="/home/kangh/lab-server/SpliceDecoder/Isoformswitch.config",
#    code_dir="/home/kangh/lab-server/SpliceDecoder/code/transcript-toolkit/",
#    species="human"
#)
##############

args, parser = parse_args(sys.argv[1:])
ref_fam = Load_data()
ref_cds = pd.read_csv(args.input+"merged_bestorf.txt",
                      sep="\t", header=None)
## Coding potential cut off
if args.species == 'human':
    cpat_cutoff = 0.364
elif args.species == 'mouse':
    cpat_cutoff = 0.44
# ref_cds = ref_cds[ref_cds[5] > cpat_cutoff]
ref_bed = pd.read_csv(args.input+"merged.bed",
                      sep="\t", header=None)

# query_list = pd.read_csv(args.query_list,
# 			 sep="\t")
# query_list.columns = ["Major","query"]

OUT = args.input+"result/"
if not os.path.exists(OUT):
    os.mkdir(OUT)

# integrity_indi = open(OUT+"Domain_integrity_indi.txt", "w")
# nmd_check = open(OUT+"NMD_check.txt", "w")
tx_gene_dict = pd.read_csv(f"{args.input}tx_gene_dict",
                           sep="\t",
                           header=None)
main_output = open(OUT+"GTF_annotation.txt", "w")
main_output.write(f"Domain\tlength(bp)\tprob_NMD\ttranscript_id\tAUG\tStop\tGene"+"\n")
# main_output.write("##Every diff is calculated by Major TX - Query TX"+"\n")
# main_output.write("Comparison"+"\t"+"Reference_transcript"+"\t"+"Query_transcript"+"\t"+"ORF"+"\t"+"AUG (Ref-Sim)"+"\t"+"Stop (Ref-Sim)"+\
#                 "\t"+"5'UTR_difference (nt)"+"\t"+"CDS_difference (nt)"+"\t"+"3'UTR_difference (nt)"+"\t"+"Domain_integrity"+"\t"+"Domain_change_rate"+\
#                 "\t"+"Length_of_reference_tx_domain"+"\t"+"Length_of_simulated_tx_domain"+"\t"+"Functional_class"+"\t"+"Probability_of_NMD"+"\n")

# %%
## Using only Top 1 ORFs
for target in ref_cds[0].unique():
    ref_start = ref_cds[ref_cds[0]==target].sort_values(by=5, ascending=False).iloc[:1,:][2]
    ref_stop = ref_cds[ref_cds[0]==target].sort_values(by=5, ascending=False).iloc[:1,:][3]
    gene_id = tx_gene_dict[tx_gene_dict[0]==target][1].values[0]
    
    def Domain_integrity(start_list, stop_list, key):
        """ Calculate all possible case to measure domain differences 

        Args:
            start_list (_type_): Start position of 3 ORFs
            stop_list (_type_): Stop position of 3 ORFs
            key (_type_): "Ref" or "Sim"

        Returns:
            _type_: total length of all possible domain
            list: total length of transcript to calculate 3'UTR diff (The end of last exon in the rel_values)
        """

        total_length = []
        bed, domain, dname = Make_query(ref_bed[ref_bed[6]==target], ref_fam[ref_fam[13]==target])
        total = pd.DataFrame()
        nmd = 0
        pre_nmd = 0
        for start, stop in zip(start_list, stop_list):  # Consider top3 ORFs of the ref_TX
            result = {}
            total_length.append(bed[-1][1])
            ####################################
            ## Updated NMD classification
            for i in range(len(bed)):
                # if bed[i][0] <= stop and \
                #    stop <= bed[i][1]:   # If the given exon contains stop codon
                #     Exon_PTC_contain = i
                #     break
                # if bed[-1][0] <= stop and \
                #     stop <= bed[-1][1]:  # Stop on a last exon
                #     nmd = nmd    # No NMD
                # elif (start - stop) > 150:   # NMD pred, it could be 50 Rule
                #     nmd = nmd    # No NMD
                # elif bed[Exon_PTC_contain][1] - bed[Exon_PTC_contain][0] > 407:   # NMD pred, it could be 50 as well
                #     nmd = nmd   # No NMD
                if (bed[-1][0] - stop) >= 55:   # NMD pred, it could be 50 as well
                    nmd += 1
                else:
                    nmd = nmd   # No NMD
            ####################################
            if start != "Loss":
                for block, domain_name in zip(domain, dname):   # Consider all domain in CDS region
                    ## Prunning domain that is out of CDS (start-stop)
                    if domain_name in result.keys():
                        did = domain_name
                    elif domain_name == " ":
                        did = did
                    else:
                        did = domain_name
                        result[did] = []
                        
                    ## Prunning steps
                    if pre_nmd == nmd:  # If it is not NMD        
                        if block[0] >= start and block[1] <= stop:
                            length = np.abs(block[1] - block[0]) + 1
                            result[did].append(length)
                            # result[did].append([block[1],block[0]])
                        elif block[0] < start and block[1] > start:
                            length = np.abs(block[1] - start) + 1
                            result[did].append(length)
                            # result[did].append([block[1],start])
                        elif block[0] < stop and block[1] > stop:
                            length = np.abs(stop - block[0]) + 1
                            result[did].append(length)
                            # result[did].append([stop,block[1]])
                    
                    else:   # If it is NMD
                        result[did].append(0.0)
                        
            result = pd.DataFrame.from_dict(result, orient="index")
            result["sum"] = result.sum(axis=1)
            result = result[["sum"]]    # make df that contains merged each same domain length (n*1) for each ORF
            
            if pre_nmd == nmd:
                result.columns = [start]    # Simplify ORF~total domain length (of each domains)
            else:
                start = str(start)+"|NMD"
                result.columns = [start]    # Simplify ORF~total domain length (of each domains)
            
            if total.shape[1] == 0:
                total = result
            else:
                total = pd.merge(total, result, 
                                 left_index=True, 
                                 right_index=True, 
                                 how="outer")

            ## pair ID, ORF, NMD_score, total_domain_length, key(Ref/Sim)
            # nmd_check.write(f"{major}_vs_{query}"+"\t"+str(start)+"\t"+str((bed[-1][0] - stop))+"\t"+str(result.sum().values[0])+"\t"+key+"\n")
            # pre_nmd = nmd   # To skip NMD form during the protein domain merge step
        
        total["Sum"] = total.sum(axis=1)
        total["key"] = key
        total["NMD"] = nmd

        return total, total_length

    ref_dom_list, ref_total_length = Domain_integrity(ref_start, ref_stop, "target")
    merged_index = ref_dom_list.index.tolist()
    comp_df = pd.DataFrame(np.vstack([ref_dom_list.to_numpy()]),
                           index=merged_index)
    
    comp_df.columns = (["ORF{}".format(orf_n+1) for orf_n in range(comp_df.shape[1]-3)]+["Sum","key","NMD"])
    comp_df = comp_df.fillna(0.0)   # Nan to zero, all loss case fill zero
    for d in comp_df.index.unique():    # All domain
        temp = comp_df[comp_df.index==d]
        if temp["Sum"].sum() == 0 and \
            len(temp["key"].unique()) == 2:
            comp_df = comp_df.drop(d)
    
    ## New for direction of Domain changes
    ## The whole_doa_direction will be chaned in line 532 depending on NMD value
    ## Annotate functional categories that have different length, it is assigned to each ORF
    # domain_only_df = comp_df[comp_df.index.str.contains(";domain;")]  # Consider only Domain or not
    domain_only_df = comp_df[~comp_df.index.str.contains(";proteome;")]
    domain_only_df = domain_only_df[['ORF1','NMD']]
    if domain_only_df.shape[0] >= 1:
        domain_only_df['NMD'] = domain_only_df['NMD'].apply(lambda x : "Yes" if x > 0 else "No")
        domain_only_df['transcript'] = target
        domain_only_df['AUG'] = ref_start.values[0]
        domain_only_df['Stop'] = ref_stop.values[0]
        domain_only_df['Gene'] = gene_id
        domain_only_df.to_csv(main_output,
                            sep="\t",
                            header=False)
            
        
main_output.close()
#nmd_check.close()
#integrity_indi.close()
# %%
