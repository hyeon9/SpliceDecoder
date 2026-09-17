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


args, parser = parse_args(sys.argv[1:])
ref_fam = Load_data()
ref_cds = pd.read_csv(args.input+"merged_bestorf.txt",
                      sep="\t", header=None)
ref_bed = pd.read_csv(args.input+"merged.bed",
                      sep="\t", header=None)
query_list = pd.read_csv(args.query_list,
                         sep="\t",
                         header=None)
query_list.columns = ["Major","query"]

OUT = args.input+"result/"
if not os.path.exists(OUT):
    os.mkdir(OUT)

integrity_indi = open(OUT+"Domain_integrity_indi.txt", "w")
nmd_check = open(OUT+"NMD_check.txt", "w")
main_output = open(OUT+"Main_output.txt", "w")

main_output.write("##Every diff is calculated by Major TX - Query TX"+"\n")
main_output.write("Comparison"+"\t"+"Reference_transcript"+"\t"+"Query_transcript"+"\t"+"ORF"+"\t"+"AUG (Ref-Sim)"+"\t"+"Stop (Ref-Sim)"+\
                "\t"+"5'UTR_difference (nt)"+"\t"+"CDS_difference (nt)"+"\t"+"3'UTR_difference (nt)"+"\t"+"Domain_integrity"+"\t"+"Domain_change_rate"+\
                "\t"+"Length_of_reference_tx_domain"+"\t"+"Length_of_simulated_tx_domain"+"\t"+"Functional_class"+"\t"+"Probability_of_NMD"+"\n")

for major in query_list["Major"].unique():
    target_query = query_list[query_list["Major"]==major]
    for query in target_query["query"].unique():
        ## Using only Top 1 ORFs
        ref_start = ref_cds[ref_cds[0]==major].sort_values(by=5, ascending=False).iloc[:1,:][2]
        ref_stop = ref_cds[ref_cds[0]==major].sort_values(by=5, ascending=False).iloc[:1,:][3]
        sim_start = ref_cds[ref_cds[0]==query].sort_values(by=5, ascending=False).iloc[:1,:][2].tolist()
        sim_stop = ref_cds[ref_cds[0]==query].sort_values(by=5, ascending=False).iloc[:1,:][3].tolist()

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
            if key == "query":
                bed, domain, dname = Make_query(ref_bed[ref_bed[6]==query], ref_fam[ref_fam[13]==query])
            else:
                bed, domain, dname = Make_query(ref_bed[ref_bed[6]==major], ref_fam[ref_fam[13]==major])
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
                            elif block[0] < start and block[1] > start:
                                length = np.abs(block[1] - start) + 1
                                result[did].append(length)
                            elif block[0] < stop and block[1] > stop:
                                length = np.abs(stop - block[0]) + 1
                                result[did].append(length)
                        
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
                nmd_check.write(f"{major}_vs_{query}"+"\t"+str(start)+"\t"+str((bed[-1][0] - stop))+"\t"+str(result.sum().values[0])+"\t"+key+"\n")
                pre_nmd = nmd   # To skip NMD form during the protein domain merge step
            
            total["Sum"] = total.sum(axis=1)
            total["key"] = key
            total["NMD"] = nmd

            return total, total_length

        ref_dom_list, ref_total_length = Domain_integrity(ref_start, ref_stop, "major")
        sim_dom_list, sim_total_length = Domain_integrity(sim_start, sim_stop, "query")
        merged_index = ref_dom_list.index.tolist() + sim_dom_list.index.tolist()
        comp_df = pd.DataFrame(np.vstack([ref_dom_list.to_numpy(),
                                          sim_dom_list.to_numpy()]),
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
        domain_only_df = comp_df[comp_df.index.str.contains(";domain;")]
        whole_doa_direction_dict = {}
        for k, orf in enumerate(comp_df.columns[:3]):   # k is for making the proper output file
            if not orf.startswith("ORF"):
                pass
            else:
                whole_doa_direction_dict[k] = []
                ref_domain_sum = domain_only_df[domain_only_df["key"]=="major"][orf].sum()
                sim_domain_sum = domain_only_df[domain_only_df["key"]=="query"][orf].sum()
                ref_other_sum = comp_df[comp_df["key"]=="major"][orf].sum()
                sim_other_sum = comp_df[comp_df["key"]=="query"][orf].sum()
                if ref_domain_sum > sim_domain_sum or \
                    ref_other_sum > sim_other_sum:
                    whole_doa_direction = "LoD"
                elif ref_domain_sum < sim_domain_sum or \
                        ref_other_sum < sim_other_sum:
                    whole_doa_direction = "GoD"
                # elif ref_other_sum != sim_other_sum:   # For non domain categories, such as motif, binding site etc
                #     whole_doa_direction = "other_regions_diff"
                else:   # No changes in any annotations
                    whole_doa_direction = "no_changes"
                whole_doa_direction_dict[k].append(whole_doa_direction)
        
        ## Calculate change ratio for each ORF
        ## It is relative change ratio
        ## 0: Not changed, 1: Fully changed
        whole_orf_diff = []
        for orf in comp_df.columns:
            if orf.startswith("ORF"):   # Sometime ORF may not be 3
                indi_orf = comp_df[[orf,"key","NMD"]]
                if (indi_orf[indi_orf["key"]=="major"][orf].sum()) == 0 and \
                   (indi_orf[indi_orf["key"]=="query"][orf].sum()) == 0:
                    indi_diff = 1.0
                
                else:   # For DOA cases
                    indi_diff = []  # For each ORF
                    for indi_dom in indi_orf.index.unique():    # Cal Change ratio for each domain
                        temp = indi_orf[indi_orf.index==indi_dom]
                        if len(temp["key"].unique()) == 2:   # If they have same domain
                            r = float(temp[temp["key"]=="major"][orf].values[0])
                            s = float(temp[temp["key"]=="query"][orf].values[0])
                            if r == 0 and s == 0:	# CDS alt
                                doa_dir = 0 # Direction of domain alteration
                                ratio = 0.0

                            else:
                                ratio = abs(s-r) / max(r,s)

                            indi_diff.append(ratio)

                            if ratio >= 0:   # Only report changed cases    , it should be > 0
                                if s > r :  # Gain of domain block
                                    doa_dir = 1
                                    integrity_indi.write(f"{major}_vs_{query}"+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                                elif s < r :    # Loss of domain block
                                    doa_dir = -1
                                    integrity_indi.write(f"{major}_vs_{query}"+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                                else:   # No changes, it won't be reported in detailed domain info file
                                    doa_dir = 0
                        
                        elif len(temp["key"].unique()) == 1:   # If the ref and sim do not have same domain
                            if temp[temp["key"]==temp["key"].unique()][orf].values[0] > 0:  # It should have real length after prunning
                                if temp["key"].unique() == "query": # If minor TX only has domain
                                    doa_dir = 1
                                    ratio = 1.0
                                else:
                                    doa_dir = -1
                                    ratio = 1.0
                                integrity_indi.write(f"{major}_vs_{query}"+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                            else:   # No changes in the domain regions, it won't be reported in detailed domain info file
                                doa_dir = 0
                                ratio = 0.0

                        elif len(temp["key"].unique()) == 0:   # If the ref and sim do not have same domain, it won't be reported in detailed domain info file
                            doa_dir = 0
                            ratio = 0.0

#                                integrity_indi.write(i+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                        indi_diff.append(ratio)
                whole_orf_diff.append(np.mean(indi_diff))    # It contains 3 different mean ratios
            
            else:
                pass
            
        ## Write output
        ref_dom_per_orf = ref_dom_list.drop(["Sum","key","NMD"], axis=1).sum()  # Sum of all functional regions (domain + motif + binding + region)
        sim_dom_per_orf = sim_dom_list.drop(["Sum","key","NMD"], axis=1).sum()
        dom_per_orf = (sim_dom_per_orf.to_numpy() / ref_dom_per_orf.to_numpy())*100 # Domain integrity
        longid = f"{major} vs {query}"
        tx = major
        occ_event = query

        for n in range(len(ref_start)):
            ref_dom_length = ref_dom_per_orf[:3].iloc[n]
            sim_dom_length = sim_dom_per_orf[:3].iloc[n]
            if sim_start[n]!="Loss":
                utr5_diff = float(ref_start.tolist()[n] - sim_start[n])
                # AA_diff = (sim_stop[n] - sim_start[n]) / \
                #           (ref_stop.tolist()[n] - ref_start.tolist()[n])    # SimTX / RefTX AA length
                AA_diff = float(ref_stop.tolist()[n] - ref_start.tolist()[n]) - \
                            float(sim_stop[n] - sim_start[n])  # RefTX AA - SimTX AA length
                utr3_diff = float(ref_total_length[n] - ref_stop.tolist()[n]) - \
                            float(sim_total_length[n] - sim_stop[n])
                dom_intig = dom_per_orf[n]
                dom_change_ratio = whole_orf_diff[n]
                if str(sim_dom_list.columns.tolist()[n]).endswith("NMD"):
                    sim_nmd = 1
                else:
                    sim_nmd = 0

                if str(ref_dom_list.columns.tolist()[n]).endswith("NMD"):
                    ref_nmd = 1
                else:
                    ref_nmd = 0

                nmd_diff = ref_nmd - sim_nmd

                
                if nmd_diff != 0:
                    final_whole_doa_direction = "NMD"
                else:   # Non-NMD, load the saved functional category information,
                    if dom_change_ratio != 0:    # There are non-zero change ratio, LoD, GoD, CDS/UTR_alt/no_changes
                        if whole_doa_direction_dict[n][0] == "no_changes":	# There are only CDS/UTR_alt/no_changes
                            if AA_diff != 0:
                                dom_change_ratio = 0	# Correction since it was defined as 1 in line 351
                                final_whole_doa_direction = "CDS_alt"
                            elif AA_diff == 0 and (utr5_diff != 0 or utr3_diff != 0):
                                dom_change_ratio = 0	# Correction since it was defined as 1 in line 351
                                final_whole_doa_direction = "UTR_alt"
                            else:
                                final_whole_doa_direction = "no_changes"
                        else:
                            final_whole_doa_direction = whole_doa_direction_dict[n][0]	# LoD or GoD
                    else:
                        final_whole_doa_direction = "no_changes"


            else:   # No matched start and stop in simTX
                utr5_diff = "Frame loss"
                AA_diff = "Frame loss"
                utr3_diff = "Frame loss"
                dom_intig = "Frame loss"
                dom_change_ratio = "Frame loss"
                final_whole_doa_direction = "Frame loss"
                nmd_diff = "Frame loss"
                
            main_output.write(longid+"\t"+
                                tx+"\t"+
                                occ_event+"\t"+
                                "pORF{}".format(n+1)+"\t"+
                                "{}-{}".format(ref_start.tolist()[n], sim_start[n])+"\t"+
                                "{}-{}".format(ref_stop.tolist()[n], sim_stop[n])+"\t"+
                                str(utr5_diff)+"\t"+
                                str(AA_diff)+"\t"+
                                str(utr3_diff)+"\t"+
                                str(dom_intig)+"\t"+
                                str(dom_change_ratio)+"\t"+
                                str(ref_dom_length)+"\t"+
                                str(sim_dom_length)+"\t"+
                                str(final_whole_doa_direction)+"\t"+
                                str(nmd_diff)+"\n")
        
main_output.close()
nmd_check.close()
integrity_indi.close()
# %%
