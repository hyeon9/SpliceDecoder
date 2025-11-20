#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import sys
import os
import argparse
from datetime import datetime

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
pd.set_option('future.no_silent_downcasting', True)
# plt.rcParams["font.family"] = "Arial"

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(description='''

Description
    #########################################################
    This script estimate functional difference between Ref_TX and Sim_TX
    This work can be focused on certain functional categoriy 
    (e.g., binding, dna_bind, motif, domain, region, and all)
    OUTPUT: _Domain_alts, _NMD_check, _Main_output
    #########################################################''',
    formatter_class=argparse.RawTextHelpFormatter)
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains rmat file', 
                        required=True,
                        type=str)
    parser.add_argument('--splicing_event', '-s', 
                        help='Splicing categories e.g., CA, RI, A3SS, A5SS, MXE', 
                        required=True,
                        type=str)
    parser.add_argument('--config_path', '-p',
                        help='The path of your config file',
                        required=True,
                        type=str)
    parser.add_argument("--code_dir", "-c", help="Path of code",
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument("--nmd_met", "-n", help="you can use other NMD finding method (default or sensitive), default: 55nt-rule", 
                        type=str,
                        default="default")
    parser.add_argument("--interesting_target", "-t", help="number of threads to use, default: all", 
                        type=str,
                        default="all")

    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args

# %%
def Load_data(key):
    sim_bed = args.input+"sim_bed/"
    cpat_DIR = args.input+"cpat/"

    ## coreID = SID|ENST|event_type
    ## Merge whole TX-splicing event pairs
    cmd = "cat "+sim_bed+"merged_"+key+ \
                "*_w_event.bed > "+ args.input+"merged_w_"+key+".bed"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    sim_tx = pd.read_csv(args.input+"merged_w_"+key+".bed",
                            sep="\t", header=None)
    sim_tx["coreID"] = sim_tx[7]+"|"+sim_tx[4]+"|"+sim_tx[6]
    
    cmd = "cat "+sim_bed+"merged_"+key+ \
                "*_wo_event.bed > "+ args.input+"merged_wo_"+key+".bed"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    ref = pd.read_csv(args.input+"merged_wo_"+key+".bed",
                            sep="\t", header=None)
    ref["coreID"] = ref[7]+"|"+ref[4]+"|"+ref[6]

    cmd = "cat "+cpat_DIR+"wo_event_bestorf_"+key+ \
            "* > "+ args.input+"wo_event_bestorf_"+key+".txt"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)
    
    cmd = "cat "+cpat_DIR+"w_event_bestorf_"+key+ \
            "* > "+ args.input+"w_event_bestorf_"+key+".txt"
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

    ## Make domain data
    cmd = f'bash {args.code_dir}/run_intersect.sh {args.config_path} {args.splicing_event}'
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

    w_pfam = pd.read_csv(args.input+"table/{}_w_Pfam.txt".format(key),
                       sep="\t", header=None)
    w_pfam["coreID"] = w_pfam[14]+"|"+w_pfam[11]+"|"+w_pfam[13]

    wo_pfam = pd.read_csv(args.input+"table/{}_wo_Pfam.txt".format(key),
                       sep="\t", header=None)
    wo_pfam["coreID"] = wo_pfam[14]+"|"+wo_pfam[11]+"|"+wo_pfam[13]

    tx_dict_file = pd.read_csv(args.input+"tx_gene_dict",
                               sep="\t", header=None)

    ################################################
    if interesting_target != "all":
        w_pfam = w_pfam[w_pfam[4]==interesting_target]
        wo_pfam = wo_pfam[wo_pfam[4]==interesting_target]
    ################################################

    return sim_tx, ref, w_pfam, wo_pfam, tx_dict_file


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
    ##############################
    strand = final_bed[5].values[0]
    pfam = pfam[pfam[5]==strand].copy()
    pfam[3] = pfam[3]+";"+pfam[4]+";"+pfam[6]
    for line in range(final_bed.shape[0]):
        if line > 0:
            if strand == "+":
                ## ADD strand specific domain assign and multiple domain case
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                intron = np.abs(int(final_bed.iloc[line-1,2] - int(final_bed.iloc[line,1])))    # (N-1)th exon end - Nth exon start
                diff.append(intron)
                new_line.append([int(final_bed.iloc[line,1])-np.sum(diff)-zero, 
                                 int(final_bed.iloc[line,2])-np.sum(diff)-zero])

                for n_domain in range(domain.shape[0]):
                    mat_start = abs(domain.iloc[n_domain][8] - domain.iloc[n_domain][1])
                    mat_end = abs(domain.iloc[n_domain][9] - domain.iloc[n_domain][2])
                    
                    if mat_start == 0 and mat_end == 0:
                        domain_start = domain.iloc[n_domain][1]
                        domain_end = domain.iloc[n_domain][2]
                    else:   # Shared region between exon and domain
                        domain_start = max(domain.iloc[n_domain][8], domain.iloc[n_domain][1])
                        domain_end = min(domain.iloc[n_domain][9], domain.iloc[n_domain][2])
                        
                    domain_line.append([int(domain_start)-np.sum(diff)-zero,
                                        int(domain_end)-np.sum(diff)-zero])
                        
                    if pre_domain == domain.iloc[n_domain,3]:
                        domain_name.append(" ")
                        pre_domain = domain.iloc[n_domain,3]
                    else:
                        domain_name.append(domain.iloc[n_domain,3])
                        pre_domain = domain.iloc[n_domain,3]

        
            else:   # Start with -2, because first exon is -1
                domain = pfam[(pfam[8]==final_bed.iloc[line,1]) &
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                intron = np.abs(np.abs(int(final_bed.iloc[line-1,1])-zero) - np.abs(int(final_bed.iloc[line,2])-zero))
                diff.append(intron)
                new_line.append([np.abs(int(final_bed.iloc[line,2])-zero)-np.sum(diff), 
                                np.abs(int(final_bed.iloc[line,1])-zero)-np.sum(diff)])    # 2*zero, represents to zero and then direction change
                
                for n_domain in range(domain.shape[0]):
                    mat_start = abs(domain.iloc[n_domain][8] - domain.iloc[n_domain][1])
                    mat_end = abs(domain.iloc[n_domain][9] - domain.iloc[n_domain][2])
                    
                    if mat_start == 0 and mat_end == 0:
                        domain_start = domain.iloc[n_domain][1]
                        domain_end = domain.iloc[n_domain][2]
                    else:   # Shared region between exon and domain
                        domain_start = max(domain.iloc[n_domain][8], domain.iloc[n_domain][1])
                        domain_end = min(domain.iloc[n_domain][9], domain.iloc[n_domain][2])
                        
                    domain_line.append([np.abs(int(domain_end)-zero)-np.sum(diff),
                                        np.abs(int(domain_start)-zero)-np.sum(diff)])

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
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                
                for n_domain in range(domain.shape[0]):
                    mat_start = abs(domain.iloc[n_domain][8] - domain.iloc[n_domain][1])
                    mat_end = abs(domain.iloc[n_domain][9] - domain.iloc[n_domain][2])
                    
                    if mat_start == 0 and mat_end == 0:
                        domain_start = domain.iloc[n_domain][1]
                        domain_end = domain.iloc[n_domain][2]
                    else:   # Shared region between exon and domain
                        domain_start = max(domain.iloc[n_domain][8], domain.iloc[n_domain][1])
                        domain_end = min(domain.iloc[n_domain][9], domain.iloc[n_domain][2])
                    
                    domain_line.append([int(domain_start)-zero,
                                        int(domain_end)-zero])
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
                              (pfam[9]==final_bed.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                
                for n_domain in range(domain.shape[0]):
                    mat_start = abs(domain.iloc[n_domain][8] - domain.iloc[n_domain][1])
                    mat_end = abs(domain.iloc[n_domain][9] - domain.iloc[n_domain][2])
                    
                    if mat_start == 0 and mat_end == 0:
                        domain_start = domain.iloc[n_domain][1]
                        domain_end = domain.iloc[n_domain][2]
                    else:   # Shared region between exon and domain
                        domain_start = max(domain.iloc[n_domain][8], domain.iloc[n_domain][1])
                        domain_end = min(domain.iloc[n_domain][9], domain.iloc[n_domain][2])
                    
                    domain_line.append([np.abs(int(domain_end)-zero),
                                        np.abs(int(domain_start)-zero)])
                        
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



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Run with -h!!")
        sys.exit(1)
    args, parser = parse_args(sys.argv[1:])
    if not args.input.endswith("/"):
        args.input = args.input+"/"

    interesting_target = args.interesting_target
    print(f"Used functional category: {interesting_target}")
        
    now = datetime.now()
    print(f"Making and Loading input from {args.splicing_event}")
    print("TIME:", now.strftime("%Y-%m-%d %H:%M:%S"))
    sim_bed, ref_bed, sim_fam, ref_fam, tx_dict_file = Load_data(args.splicing_event)
    ref_cds = pd.read_csv(args.input+"wo_event_bestorf_{}.txt".format(args.splicing_event),
                          sep="\t", header=None)
    sim_cds = pd.read_csv(args.input+"w_event_bestorf_{}.txt".format(args.splicing_event),
                          sep="\t", header=None)

    # ref_cds["coreID"] = ref_cds[0]+"|"+ref_cds[1].str.split("_").str[0]+"|"+ref_cds[2]
    # sim_cds["coreID"] = sim_cds[0]+"|"+sim_cds[1].str.split("_").str[0]+"|"+sim_cds[2]
    ref_cds["coreID"] = ref_cds[0]+"|"+ref_cds[1].str.replace("_wo_event","")+"|"+ref_cds[2]
    sim_cds["coreID"] = sim_cds[0]+"|"+sim_cds[1].str.replace("_w_event","")+"|"+sim_cds[2]
    ## coreID = SID|ENST|event_type

    import os
    OUT = args.input+"result/"
    if not os.path.exists(OUT):
        os.mkdir(OUT)

    integrity_indi = open(OUT+"{}_{}_Domain_alts.tsv".format(interesting_target, args.splicing_event), "w")
    nmd_check = open(OUT+"{}_{}_NMD.tsv".format(interesting_target, args.splicing_event), "w")
    main_output = open(OUT+"{}_{}_Main_table.tsv".format(interesting_target,args.splicing_event), "w")
    
    integrity_indi.write("Pair_info\tpORF\taltered_domain\tDomain_change_ratio\tChange_direction\n")
    nmd_check.write("LongID\tAUG\tpORF\tdistance(last_exon_junction-stop)\ttotal_domain_length\tkey(Ref/Sim)\tNMD_possibility\tcontain_PTC\n")
    main_output.write("##Every diff is calculated by Ref TX - Sim TX"+"\n")
    # main_output.write("LongID"+"\t"+"Target_TX"+"\t"+"occurred_event"+"\t"+"ORF_priority"+"\t"+"Start"+"\t"+"Stop"+\
    #                 "\t"+"5'UTR"+"\t"+"dAA"+"\t"+"3'UTR"+"\t"+"Domain_integrity"+"\t"+"Domain_change_ratio"+\
    #                 "\t"+"Ref_domain_length"+"\t"+"Sim_domain_length"+"\t"+"Functional_class"+"\t"+"pNMD"+"\n")
    main_output.write("LongID"+"\t"+"Gene symbol"+"\t"+"Reference_transcript"+"\t"+"Simulated_event"+"\t"+"ORF"+"\t"+"AUG (Ref-Sim)"+"\t"+"Stop (Ref-Sim)"+\
                    "\t"+"5'UTR_difference (nt)"+"\t"+"CDS_difference (nt)"+"\t"+"3'UTR_difference (nt)"+"\t"+"Domain_integrity"+"\t"+"Domain_change_rate"+\
                    "\t"+"Length_of_reference_tx_domain"+"\t"+"Length_of_simulated_tx_domain"+"\t"+"Functional_class"+"\t"+"Probability_of_NMD"+"\n")
    
    now = datetime.now()
    print(f"Performing a {args.splicing_event} simulation...")
    print("TIME:", now.strftime("%Y-%m-%d %H:%M:%S"))
    for i in ref_bed["coreID"].unique():    ## Specific sim result
        if i in ref_fam["coreID"].tolist(): ## Only considered protein coding tx
            def Trans_code(orf_list):
                sub_ref_bed = ref_bed[ref_bed["coreID"]==i].copy()
                sub_ref_bed["exon_length"] = np.abs(sub_ref_bed[2]) - np.abs(sub_ref_bed[1])    # Dose not need +1, when make this bed flank start-1 (in 01-4*.py)
                sub_sim_bed = sim_bed[sim_bed["coreID"]==i].copy()
                sub_sim_bed["exon_length"] = np.abs(sub_sim_bed[2]) - np.abs(sub_sim_bed[1])
                
                def Cal():
                    sim_orf_list = []
                    for ref_start_pos in orf_list:  # Consider all ORFs
                        ## Decoder
                        position = np.array([])
                        for line in range(sub_ref_bed.shape[0]):    # Per exon
                            position = np.append(position, sub_ref_bed.iloc[line,-1])    # Cumulate exon length
                            if ref_start_pos <= (position.sum()):   # For find exon that contains start codon
                                if sub_ref_bed[5].values[0] == "-":
                                    genomic_pos = (sub_ref_bed.iloc[line,2]-(ref_start_pos - position[:-1].sum()))
                                    break
                                else:
                                    genomic_pos = (sub_ref_bed.iloc[line,1]+(ref_start_pos - position[:-1].sum()))
                                    ## The start of the exon (Genomic coord) + (end of the exon (Genomic coord))
                                    ## Why does minus position[:-1].sum() what's the meaning???
                                    break
                        
                        ## Encoder
                        position = np.array([])
                        start_codon = "Loss"
                        for line in range(sub_sim_bed.shape[0]):
                            position = np.append(position, sub_sim_bed.iloc[line,-1])
                            if genomic_pos >= sub_sim_bed.iloc[line,1] and \
                            genomic_pos <= sub_sim_bed.iloc[line,2]:    # Check if the Ref-TX start/stop site is on the Sim-TX exon
                                if sub_sim_bed[5].values[0] == "-":
                                        start_codon = (sub_sim_bed.iloc[line,2]-genomic_pos)+position[:-1].sum()
                                        break
                                else:
                                        start_codon = (genomic_pos-sub_sim_bed.iloc[line,1])+position[:-1].sum()
                                        break
                        
                        sim_orf_list.append(start_codon)

                    return sim_orf_list
                
                sim_orf_list = Cal()

                return sim_orf_list

            ref_start = ref_cds[ref_cds["coreID"]==i].sort_values(by=7, ascending=False).iloc[:3,:][4]  # Using only top3 ORFs
            ref_stop = ref_cds[ref_cds["coreID"]==i].sort_values(by=7, ascending=False).iloc[:3,:][5]
            if ref_cds[ref_cds["coreID"]==i].shape[0] == 0:
                print(f"Check your transcript ID {i}")

            try:
                temp_sim_start = Trans_code(ref_start)   # Find a same position of ref ORF in sim tx
                temp_sim_stop = Trans_code(ref_stop)   # Find a same position of ref ORF in sim tx
            except:
                now = datetime.now()
                print(i)
                print("TIME:", now.strftime("%Y-%m-%d %H:%M:%S"))
            
            sim_start = []
            sim_stop = []
            check = []
            sub_sim_cds = sim_cds[sim_cds["coreID"]==i] # Predicted ORF in sim TX
            for rank,sim_st in enumerate(temp_sim_start):
                if sim_st in sub_sim_cds[4].tolist():   # If appORF exist in sim_TX ORF data
                    sim_start.append(sim_st)
                    sim_stop.append(sub_sim_cds[sub_sim_cds[4]==sim_st][5].values[0])

                elif sim_st == 'Loss':  # There are not matched start codon on the genomic position due to the AS 
                    if temp_sim_stop[rank] in sub_sim_cds[5].tolist():   # different start, same end (partial frame)
                        sim_start.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][4].values[0])
                        sim_stop.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][5].values[0])
                    else:   # No matched start and stop in Sim TX
                        sim_start.append("Loss")    # Frame loss
                        sim_stop.append(500)    # It is considered as NMD

                else:   # Despite there are matched genomic postion, they lose their CP
                    if temp_sim_stop[rank] in sub_sim_cds[5].tolist():   # different start, same end (partial frame)
                        if not i.startswith("MXE"): # Check length of DS event
                            l_DS = np.abs(int(i.split(";")[6]) - int(i.split(";")[7]))
                            if (l_DS+1)%3 == 0: ## Updated
                                sim_start.append(sim_st)
                                sim_stop.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][5].values[0])
                            else:
                                sim_start.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][4].values[0])
                                sim_stop.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][5].values[0])
                        else:
                            l_DS = np.abs(int(i.split(";")[6]) - int(i.split(";")[7]))
                            l_DS2 = np.abs(int(i.split(";")[8]) - int(i.split(";")[9]))
                            if (l_DS+1)%3 == 0 and (l_DS2+1)%3 == 0:    ## Updated
                                sim_start.append(sim_st)
                                sim_stop.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][5].values[0])
                            else:
                                sim_start.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][4].values[0])
                                sim_stop.append(sub_sim_cds[sub_sim_cds[5]==temp_sim_stop[rank]][5].values[0])
                    else:   # No matched start and stop in Sim TX
                        sim_start.append("Loss")    # Frame loss
                        sim_stop.append(500)    # It is considered as NMD
    

            def Domain_integrity(start_list, stop_list, key):
                """ Consider all reading frame to measure domain differences 

                Args:
                    start_list (_type_): Start position of 3 ORFs
                    stop_list (_type_): Stop position of 3 ORFs
                    key (_type_): "Ref" or "Sim"

                Returns:
                    _type_: total length of all possible domain
                    list: total length of transcript to calculate 3'UTR diff (The end of last exon in the rel_values)
                """
                
                total_length = []
                if key == "Sim":
                    bed, domain, dname = Make_query(sim_bed[sim_bed["coreID"]==i], sim_fam[sim_fam["coreID"]==i])
                else:
                    bed, domain, dname = Make_query(ref_bed[ref_bed["coreID"]==i], ref_fam[ref_fam["coreID"]==i])
                total = pd.DataFrame()
                nmd = 0
                pre_nmd = 0
                orf_num = 0

                for start, stop in zip(start_list, stop_list):  # Consider top3 ORFs of the ref_TX
                    ptc_contain = 0
                    orf_num += 1
                    result = {}
                    total_length.append(bed[-1][1])
                    nmd_annot = "Loss"
                    if start != "Loss":
                        start = int(start)
                        stop = int(stop)
                        
                        if args.nmd_met == "default":   # It only consider 55 rule
                            if (bed[-1][0] - stop) >= 55:   # NMD pred, it could be 50 as well
                                nmd += 1
                                # nmd_annot = "NMD (50-55rule)"
                                nmd_annot = "HIGH (55nt)"
                                ptc_contain += 1
                            else:
                                nmd = nmd   # No NMD
                                # nmd_annot = "No NMD"
                                nmd_annot = "No"
                                
                        else:   # Applying Nat Genet paper criteria (sensitive)
                            for exon_n in range(len(bed)):   # Find the exon that contains stop codon
                                if bed[exon_n][0] <= stop and \
                                    stop <= bed[exon_n][1]:   # If the given exon contains stop codon
                                    Exon_PTC_contain = exon_n
                                    break
                            
                            if stop < bed[-1][0]: # stop codon locates not in last exon == PTC
                                ptc_contain += 1
                                if abs(start - stop) < 150: # Start-proximal evade NMD
                                    nmd = nmd
                                    # nmd_annot = "No NMD (Start-proximal)"
                                    nmd_annot = "INTERMEDIATE (Start-proximal)"
                                elif abs(bed[Exon_PTC_contain][1] - bed[Exon_PTC_contain][0]) > 407:   # Long exon evade NMD
                                    nmd = nmd
                                    # nmd_annot = "No NMD (Long-exon)"
                                    nmd_annot = "INTERMEDIATE (Long-exon)"
                                elif (bed[-1][0] - stop) >= 55:   # NMD pred, it could be 50 as well
                                    nmd += 1
                                    # nmd_annot = "NMD (50-55rule)"
                                    nmd_annot = "HIGH (55nt)"
                                else:
                                    nmd = nmd   # No NMD   
                                    # nmd_annot = "No NMD"
                                    nmd_annot = "LOW (less 55nt)"
                            else:   # PTC in last exon
                                ptc_contain = 0
                                nmd = nmd   # No NMD
                                # nmd_annot = "No NMD"
                                nmd_annot = "No"
                     
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

                    ## NMD TABLE: pair ID, ORF, NMD_score, total_domain_length, key(Ref/Sim), and NMD_type
                    if ptc_contain > 0:
                        ptc_class = 'Y'
                    else:
                        ptc_class = 'N'
                    nmd_check.write(i+"\t"+str(start)+"\t"+f'ORF{orf_num}'+"\t"+str((bed[-1][0] - stop))+"\t"+str(result.sum().values[0])+"\t"+key+"\t"+nmd_annot+"\t"+str(ptc_class)+"\n")
                    pre_nmd = nmd   # To skip NMD form during the protein domain merge step
                
                total["Sum"] = total.sum(axis=1)
                total["key"] = key
                total["NMD"] = nmd

                return total, total_length


            ref_dom_list, ref_total_length = Domain_integrity(ref_start, ref_stop, "Ref")
            sim_dom_list, sim_total_length = Domain_integrity(sim_start, sim_stop, "Sim")
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
                    ref_domain_sum = domain_only_df[domain_only_df["key"]=="Ref"][orf].sum() 
                    sim_domain_sum = domain_only_df[domain_only_df["key"]!="Ref"][orf].sum()
                    ref_other_sum = comp_df[comp_df["key"]=="Ref"][orf].sum() 
                    sim_other_sum = comp_df[comp_df["key"]!="Ref"][orf].sum()
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
            ## Each ORF has functional category that has different length and their different length information
            whole_orf_diff = []
            for orf in comp_df.columns:
                if orf.startswith("ORF"):   # Sometime ORF may not be 3
                    indi_orf = comp_df[[orf,"key","NMD"]]   # Consider all categories
                    if (indi_orf[indi_orf["key"]=="Ref"][orf].sum()) == 0 and \
                       (indi_orf[indi_orf["key"]=="Sim"][orf].sum()) == 0:  # Non functional transcripts
                        indi_diff = 0.0
                    
                    elif (indi_orf[indi_orf["key"]=="Ref"][orf].sum()) == 0 or \
                         (indi_orf[indi_orf["key"]=="Sim"][orf].sum()) == 0:  # For NMD related cases
                        indi_diff = 1.0
                    
                    else:   # For DOA cases >> integrity_indi only has DOA cases, it does not have NMD cases
                        indi_diff = []  # For each ORF
                        for indi_dom in indi_orf.index.unique():    # Cal Change ratio for each domain block
                            temp = indi_orf[indi_orf.index==indi_dom]   # Each domain block
                            
                            if len(temp["key"].unique()) == 2:   # If they have same domain
                                r = float(temp[temp["key"]=="Ref"][orf].values[0])  # Whole length of certain domain of refTX
                                s = float(temp[temp["key"]=="Sim"][orf].values[0])  # Whole length of certain domain of simTX
   
                                if r == 0 and s == 0:
                                    doa_dir = 0 # Direction of domain alteration
                                    ratio = 0.0
                                
                                else:
                                    ratio = abs(s-r) / max(r,s)

                                indi_diff.append(ratio)

                                if ratio >= 0:   # Only report changed cases    , it should be > 0
                                    if s > r :  # Gain of domain block
                                        doa_dir = 1
                                    elif s < r :    # Loss of domain block
                                        doa_dir = -1
                                    else:   # No changes    IF main if is > 0, this part should be removed
                                        doa_dir = 0
                                    integrity_indi.write(i+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                                    
                            elif len(temp["key"].unique()) == 1:   # If the ref and sim do not have same domain
                                if temp[temp["key"]==temp["key"].unique()][orf].values[0] > 0:  # It should have real length after prunning
                                    if temp["key"].unique() == "Sim":
                                        doa_dir = 1
                                        ratio = 1.0
                                    else:
                                        doa_dir = -1
                                        ratio = 1.0
                                    integrity_indi.write(i+"\t"+str(orf)+"\t"+indi_dom+"\t"+str(ratio)+"\t"+str(doa_dir)+"\n")  # Save individual informaion for further analysis
                                else:
                                    doa_dir = 0
                                    ratio = 0.0
                                    
                            elif len(temp["key"].unique()) == 0:   # If the ref and sim do not have same domain
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
            longid = i.split("|")[0]
            tx = i.split("|")[1]
            occ_event = i.split("|")[2]

            for n in range(len(ref_start)):
                ref_dom_length = ref_dom_per_orf[:3].iloc[n]
                sim_dom_length = sim_dom_per_orf[:3].iloc[n]
                if sim_start[n]!="Loss":
                    utr5_diff = float(ref_start.tolist()[n] - sim_start[n])
                    AA_diff = float(ref_stop.tolist()[n] - ref_start.tolist()[n]) - \
                              float(sim_stop[n] - sim_start[n])  # RefTX AA - SimTX AA length 
                    utr3_diff = float(ref_total_length[n] - ref_stop.tolist()[n]) - \
                                float(sim_total_length[n] - sim_stop[n])
                    dom_intig = dom_per_orf[n]
                    dom_change_ratio = whole_orf_diff[n]
                    if "|NMD" in str(sim_dom_list.columns.tolist()[n]):
                        sim_nmd = 1
                    else:
                        sim_nmd = 0
                    
                    if "|NMD" in str(ref_dom_list.columns.tolist()[n]):
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

                gene_symbol = str(tx_dict_file[tx_dict_file[0]==tx][2].unique()[0])
                main_output.write(longid+"\t"+
                                  gene_symbol+"\t"+
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

now = datetime.now()
print("TIME:", now.strftime("%Y-%m-%d %H:%M:%S"))
print(f'{args.splicing_event} is DONE')
print(f'Your output folder: {OUT}')
# %%
