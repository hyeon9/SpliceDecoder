#!/usr/bin/env python3
# %%
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import subprocess
import argparse
import os
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Arial"
plt.rcParams.update({
'axes.titlesize': 15,     # 제목 글꼴 크기
'axes.labelsize': 14,     # x, y축 라벨 글꼴 크기
'xtick.labelsize': 12,    # x축 틱 라벨 글꼴 크기
'ytick.labelsize': 12,    # y축 틱 라벨 글꼴 크기
'legend.fontsize': 12,    # 범례 글꼴 크기
'figure.titlesize': 15    # figure 제목 글꼴 크기
})

def parse_args(cmd_args=None, namespace=None):
    parser=argparse.ArgumentParser(
        description='''This script discoveries match isoforms with given differential splicing events''')
    ## Positional
    parser.add_argument('--input', '-i', 
                        help='Input directory that contains simulated results', 
                        required=True,
                        type=str)
    parser.add_argument('--cano_tx', '-g', 
                        help='Canonical_transcript ID',
                        required=True,
                        type=str)
    
    ## Optional
    parser.add_argument('--remove_info', '-ri', nargs='*',
                        help='DO NOT DEPIC these functional category e.g., coiled chain', 
                        type=str,
                        default=[])                   
    
    args = parser.parse_args(cmd_args, namespace)
    
    return args, parse_args
    

args, parser = parse_args(sys.argv[1:])
args.input = args.input+"/"

################################################
################################################ TEST RUN
# if __name__ == "__main__":
#     test_args = [
#         "--input", "/home/kangh/lab-server/KAIST/Tx_seq_compare",
#         "--cano_tx", "ENST00000286186.11",
#         "--query_list", "/home/kangh/lab-server/KAIST/Tx_seq_compare/Seq_comparison_template.tsv",
#         "-ri", "region",
#     ]
#     args, parser = parse_args(test_args)
#     args.input = args.input + "/"

#     print(args)
################################################
# %%
def Run(key):
    
    def Load_data(key, remove_info):
        """ Load information data
            It could be smaller (perhaps using cmd?)
        Args:
            key (_type_): Cano-TX

        Returns:
            _type_: _description_
        """
        
        query_list = pd.read_csv(args.input+"query_list.txt",
                                 sep="\t")
        query_list.columns = ["Major","query"]
        query_list = query_list[query_list["Major"]==key]   # Fix bug
        merged_bed = pd.read_csv(args.input+"merged.bed",
                                 sep="\t",
                                 header=None)
        ref_tx = merged_bed[merged_bed[6].isin(query_list["Major"].tolist())]
        sim_tx = merged_bed[merged_bed[6].isin(query_list["query"].tolist())]
        query = pd.read_csv(args.input+"result/Main_output.txt",
                            sep="\t", skiprows=1)
        
        query = query[query["5'UTR_difference (nt)"]!="Frame loss"] # Filter out frame loss, because they donot have domain
        query = query[query["ORF"]=="pORF1"]
        query["ORF"] = query["ORF"].apply(lambda x : 3 if x =="pORF1" else (2 if x == "pORF2" else 1))
        query = query[["Comparison","Query_transcript","AUG (Ref-Sim)","Stop (Ref-Sim)","CDS_difference (nt)","Reference_transcript","ORF","Domain_change_rate","Functional_class","Probability_of_NMD"]]
        query["simStop"] = query["Stop (Ref-Sim)"].str.split("-").str[1]
        query["Stop (Ref-Sim)"] = query["Stop (Ref-Sim)"].str.split("-").str[0]
        query["simStart"] = query["AUG (Ref-Sim)"].str.split("-").str[1]
        query["AUG (Ref-Sim)"] = query["AUG (Ref-Sim)"].str.split("-").str[0]
        query.columns = ["ID","query_tx","Start","STOP","dAA","Major_tx","ORF","delta L","DOA_types","pNMD","simSTOP","simStart"]
        query["comp_pair"] = query["ID"]
        query[["simSTOP","STOP","Start","simStart"]] = query[["simSTOP","STOP","Start","simStart"]].astype(float)
        query = query[query['Major_tx']==args.cano_tx].copy()

        pfam = pd.read_csv(args.input+"/table/Pfam.txt",
                                 sep="\t",
                                 header=None)
        
        wo_pfam = pfam[pfam[13].isin(query_list["Major"].tolist())]
        w_pfam = pfam[pfam[13].isin(query_list["query"].tolist())]
    
        ## Filter out required info
        w_pfam = w_pfam[~w_pfam[4].isin(remove_info)]
        wo_pfam = wo_pfam[~wo_pfam[4].isin(remove_info)]

        w_pfam[3] = w_pfam[3] +" ("+w_pfam[4]+")"
        wo_pfam[3] = wo_pfam[3] +" ("+wo_pfam[4]+")"

        return sim_tx, ref_tx, query, w_pfam, wo_pfam


    def Comp_plot(canvas, exonpos, exon_color, *pos_info):
        """ Draw comparison plot

        Args:
            canvas (_type_): subplots axs
            exonpos (list): exon position
            exon_color (_type_): exon block color
            pos_info (list): start, stop_codon, and different exon
        Returns:
            _type_: Figure
        """
        sns.despine(left=True, top=True, bottom=True)
        canvas.tick_params(left=False, bottom=False)
        canvas.set_xticks([])
        canvas.set_yticks([])
        bgcolor="#FFFFFF"  # background
        if len(exonpos) > 0:  # If there are available exon coordinates
            ylims = {'exon_max':2.0, 'exon_min':.5}   # The height of exon box
            canvas.set_facecolor(bgcolor)
            
            def set_limits():
                ylims['bar_min'] = ylims['exon_max']+0.2
                ylims['bar_max'] = ylims['bar_min']+(ylims['exon_max']-ylims['exon_min'])/5.0

            def draw_exon(span, canvas, *domain):
                """ Draw exon

                Args:
                    span (_type_): [exon start, exon end]
                    canvas (_type_): axs
                    *domain (list): Start, Stop codon, and changed postion, len(domain) >2 means they have target exon
                """
                
                ## Mark Start and Stop codon
                if len(domain) >= 2:
                    scale_factor_max = 0.7 
                    scale_factor_min = 2
                    if span[0] < domain[0] and span[1] < domain[0]: # It is totally 5' UTR 
                        canvas.fill_between(span, 
                                            ylims['exon_min']*scale_factor_min,
                                            ylims['exon_max']*scale_factor_max,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                    
                    elif span[0] < domain[0] and span[1] >= domain[0]: # It is partially 5' UTR
                        canvas.fill_between([span[0],domain[0]], 
                                            ylims['exon_min']*scale_factor_min, 
                                            ylims['exon_max']*scale_factor_max,
                                            edgecolor="none", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        canvas.fill_between([domain[0],span[1]], 
                                            ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="#000000", facecolor=exon_color,
                                            zorder=4, linewidth=0.5)
                        ## Draw outlines
                        canvas.plot([span[0], domain[0]], [ylims['exon_min']*scale_factor_min, ylims['exon_min']*scale_factor_min], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], domain[0]], [ylims['exon_max']*scale_factor_max, ylims['exon_max']*scale_factor_max], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[0], span[1]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[0], span[1]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], span[0]], [ylims['exon_min']*scale_factor_min, ylims['exon_max']*scale_factor_max], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[1], span[1]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        
                    elif span[0] < domain[1] and span[1] > domain[1]:   # It is partially 3' UTR
                        canvas.fill_between([span[0],domain[1]], 
                                            ylims['exon_min'], 
                                            ylims['exon_max'],
                                            edgecolor="#000000", facecolor=exon_color,
                                            zorder=4, linewidth=0.5)
                        canvas.fill_between([domain[1],span[1]], 
                                            ylims['exon_min']*scale_factor_min, 
                                            ylims['exon_max']*scale_factor_max,
                                            edgecolor="none", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        ## Draw outlines
                        canvas.plot([span[0], domain[1]], [ylims['exon_min'], ylims['exon_min']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], domain[1]], [ylims['exon_max'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[1], span[1]], [ylims['exon_min']*scale_factor_min, ylims['exon_min']*scale_factor_min], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([domain[1], span[1]], [ylims['exon_max']*scale_factor_max, ylims['exon_max']*scale_factor_max], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[0], span[0]], [ylims['exon_min'], ylims['exon_max']], color="#000000", linewidth=0.5, zorder=4)
                        canvas.plot([span[1], span[1]], [ylims['exon_min']*scale_factor_min, ylims['exon_max']*scale_factor_max], color="#000000", linewidth=0.5, zorder=4)

                    elif span[0] > domain[1]:   # It is totally 3' UTR
                        canvas.fill_between(span, 
                                            ylims['exon_min']*scale_factor_min, 
                                            ylims['exon_max']*scale_factor_max,
                                            edgecolor="#000000", facecolor="#D8D8D8",
                                            zorder=3, linewidth=0.5)
                        
                    else:   # CDS
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor=exon_color,
                                            zorder=3, linewidth=0.5)
                    
                    ########################################################################
                
                ## Draw domain block
                elif len(domain) == 1:
                    canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor="#FFFFFF", facecolor=exon_color,
                                        zorder=3, linewidth=0.3)
                    canvas.text(0.9,0.26,"{}".format(domain[0]),
                                fontsize=8, zorder=4, 
                                transform=canvas.transAxes)   # Add domain name
                else:
                    canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                        edgecolor="#FFFFFF", facecolor=exon_color,
                                        zorder=3, linewidth=0.3)

                return True


            def draw_plot(canvas):
                set_limits()
                ## draw_exon line : 149
                domain_num = 0  # Set start point
                for i in range(len(exonpos)):   # For each exon block
                    ## Drawing the protein domain
                    if len(pos_info) == 1:  # If Sim,Ref-TX has CDS differences
                        if domain_num == 0:
                            draw_exon(exonpos[i], canvas, pos_info[0][domain_num])
                        else:
                            draw_exon(exonpos[i], canvas)
                        domain_num += 1

                    ## Drawing the unchanged exon structure
                    elif len(pos_info) > 1 and len(pos_info) <= 2:
                        draw_exon(exonpos[i], canvas, *pos_info)

                canvas.fill_between([exonpos[0][0], exonpos[-1][1]],
                                    ylims['bar_min'], ylims['bar_max'],
                                    edgecolor=bgcolor, facecolor=bgcolor)

            draw_plot(canvas)
            

    def Make_query(final, pfam):
        """ Collect genomic coordinate information from query event

        Args:
            final (dataframe): 01-4 output (merged_bed)

        Returns:
            new_line : Input of the drawing figure function
        """
        diff = []
        new_line = []
        domain_line = []
        domain_name = []

        for i in range(final.shape[0]):
            exon_start = int(final.iloc[i, 1])
            exon_end = int(final.iloc[i, 2])
            exon_len = abs(exon_start - exon_end)
            diff.append(exon_len)

            if i == 0:
                new_line.append([0, exon_len])
            else:
                new_line.append([new_line[i-1][1], new_line[i-1][1] + exon_len])

            domain = pfam[(pfam[8] == exon_start) & (pfam[9] == exon_end)][[0, 1, 2, 3, 5]].drop_duplicates()

            if domain.shape[0] > 0:
                n_domain_start = int(domain.iloc[0, 1])
                n_domain_end = int(domain.iloc[0, 2])
                domain_label = domain.iloc[0, 3]

                exon_cum_start = new_line[i][0]
                
                def make_domain(d_start, d_end, e_start, e_end, strand):
                    mat_start = abs(d_start - e_start)
                    mat_end = abs(d_end - e_end)
                    
                    if mat_start == 0 and mat_end == 0:
                        domain_start = d_start
                        domain_end = d_end
                    else:   # Shared region between exon and domain
                        domain_start = max(d_start, e_start)
                        domain_end = min(d_end, e_end)
                    
                    if strand == "+":
                        if domain_start != e_start:
                            rel_start = abs(domain_start - e_start) + exon_cum_start
                        else:
                            rel_start = exon_cum_start
                    
                    else:    
                        if domain_end != e_end:
                            rel_start = abs(domain_end - e_end) + exon_cum_start
                        else:
                            rel_start = exon_cum_start
                        
                    return rel_start, domain_start, domain_end

                rel_start, domain_start, domain_end = make_domain(n_domain_start, n_domain_end, exon_start, exon_end, domain[5].unique()[0])
                rel_end = rel_start + abs(domain_start - domain_end)
                domain_line.append([rel_start, rel_end])

                if domain_label not in domain_name:
                    domain_name.append(domain_label)
                else:
                    domain_name.append(" ")
        
        return new_line, domain_line, domain_name


    sim_tx, ref, query, w_pfam, wo_pfam = Load_data(key, args.remove_info)
    test_list = query.iloc[:,0].unique()
    comp_pair = test_list[0].split("|")[0]
    test_list = query[query["Major_tx"].str.contains(args.cano_tx)]["query_tx"].unique()   # Using Canonical TX ID to extract target results
    
    ## Remove Biding site in figures, it spends lots of space
    ## Capture target domains which are contained in target transcript
    # w_pfam = w_pfam[(w_pfam[4]!="binding")]
    # wo_pfam = wo_pfam[(wo_pfam[4]!="binding")]
    n_domain = wo_pfam[3].unique()    # If the given Cano-tx has at lesat one domain, the figure will be generated
    
    OUT = args.input+"figure/"
    if not os.path.exists(OUT):
        os.mkdir(OUT)
        
    for id in test_list:    # Certain splicing event with one of the possible ORF
        LID = id.split("|")[0]
        sub = query[query["query_tx"]==id]
        query_pfam = w_pfam[w_pfam[13]==id]

        ######## Functional Annotation
        doa_types = str(sub["DOA_types"].values[0])
        if sub["pNMD"].values[0] == "-1":
            doa_types = "NMD"
        elif sub["pNMD"].values[0] == "1":
            doa_types = "PTC removal"
        ########

        deltaL = float(sub["delta L"].values[0])
        sub_sim_tx = sim_tx[(sim_tx[6]==sub["query_tx"].values[0])].copy()
        sub_ref = ref[(ref[6]==sub["Major_tx"].values[0])].copy()
        
        def Draw_splicemap(fig_input,ref,sim_tx,fig_on):
            if fig_input.shape[0] == 0: # If input is None
                return
            
            else:
                n_row = 0   # For TX
                diff_dom = 0    # For domain
                fig = plt.figure(figsize=(8,(len(n_domain)/3)+1.4))   # Main figure template
                # fig = plt.figure(figsize=(8,3))
                
                for num, kinds_domain in enumerate(range(len(n_domain))):   # Make Figure and dataframe "one domain at once"
                    ## Prepare figure params
                    gene = fig_input["ID"].values[0]
                    tx = fig_input["Major_tx"].values[0]
                    eventid = fig_input["query_tx"].values[0]

                    ## Make an input (bed) for visualization
                    ref_bed, ref_domain, ref_dname = Make_query(ref, wo_pfam[wo_pfam[3]==n_domain[kinds_domain]])
                    sim_bed, sim_domain, sim_dname = Make_query(sim_tx, query_pfam[query_pfam[3]==n_domain[kinds_domain]])

                    ################################
                    ## Make a prunned bed
                    ## Make a dataframe for subsequent analysis
                    def Prunning(domain_bed, dname, st, end):
                        ## Add truncated domain
                        ## Motif could be located on UTR regions
                        l_domain = []   # length
                        p_domain_bed = []   # domain position
                        p_dname = []    # domain name
                        for index in range(len(domain_bed)):
                            ## Prunning
                            if domain_bed[index][0] >= st and \
                                domain_bed[index][1] <= end:
                                l_domain.append((domain_bed[index][1] - domain_bed[index][0]))
                                p_domain_bed.append(domain_bed[index])
                                p_dname.append(dname[index])

                            elif domain_bed[index][0] < end and \
                                domain_bed[index][1] > end: # Partial domain
                                l_domain.append((end - domain_bed[index][0]))
                                p_domain_bed.append([domain_bed[index][0],end])
                                p_dname.append(dname[index])

                            elif domain_bed[index][0] < st and \
                                domain_bed[index][1] > st: # Partial domain
                                l_domain.append((domain_bed[index][1] - st))
                                p_domain_bed.append([st,domain_bed[index][1]])
                                p_dname.append(dname[index])

                        return p_domain_bed, p_dname, l_domain

                    p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                       fig_input["Start"].values[0], 
                                                                       fig_input["STOP"].values[0])
                    p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                       fig_input["simStart"].values[0], 
                                                                       fig_input["simSTOP"].values[0])

                    ###################################
                    if (len(p_ref_domain) > 0 or len(p_sim_domain) > 0) and fig_on == "on":
                        if n_row == 0:
                            ## Main complex figure
                            ref_axs = plt.subplot(2+(len(n_domain)*2), 1, 1)
                                ## If sim splicing is inclusion form, ref tx should have same tx structure
                                ## Draw Ref-TX
                            Comp_plot(ref_axs, ref_bed, "#FF7878", 
                                      fig_input["Start"].values[0], 
                                      fig_input["STOP"].values[0])
                            plt.text(-0.08,0.2,"{}".format("Cano-Tx"), 
                                    transform=ref_axs.transAxes,
                                    fontsize=10)
                            
                            value_dict = {"PTC removal":"NMD",
                                        "NMD":"Intact",
                                        "LoD":"GoD",
                                        "GoD":"LoD",
                                        "UTR_alt":"",
                                        "CDS_alt":"",
                                        "no_changes":"No Changes"} # BUG fix 25.01.27
                            ref_axs.text(.96,.2,f"{value_dict[doa_types]}",
                                        fontsize=10, zorder=5,
                                        transform=ref_axs.transAxes)
                            
                            ## Draw Sim-TX
                            sim_axs = plt.subplot(2+(len(n_domain)*2), 1, 2, sharex=ref_axs)
                            Comp_plot(sim_axs, sim_bed, "#73ABFF", 
                                        fig_input["simStart"].values[0], 
                                        fig_input["simSTOP"].values[0])
                                                                
                            plt.text(-0.08,0.2,"{}".format("Query-Tx"), 
                                    transform=sim_axs.transAxes,
                                    fontsize=10)
                            sim_axs.text(.96,.2,f"{doa_types}",
                                        fontsize=10, zorder=5,
                                        transform=sim_axs.transAxes)
                            
                            n_row += 1  # For next domain, tx drawing function will be skipped

                        ## Draw protein domain
                        ## Last exon-exon junction - stop codon, it may be less than 55 it will be concomp_pairered as NMD
                        ref_nmd = ref_bed[-1][0] - fig_input["STOP"].values[0]
                        sim_nmd = sim_bed[-1][0] - fig_input["simSTOP"].values[0]
                        
                        def Domain_plot(diff_dom):
                            axs2 = plt.subplot(2+(len(n_domain)*2), 1, 3+diff_dom, sharex=ref_axs)
                            ## Check here, domain block concordance
                            ## It should be concomp_pairered skipped or included exon's length. It is key!!!!
                            lod = []
                            god = []
                            
                            ## In case of one of transcript does not have domain
                            if len(p_ref_domain) == 0 and len(p_sim_domain) != 0:
                                for i in range(len(p_sim_domain)):
                                    god.append([p_sim_domain[i][0], p_sim_domain[i][1]])
                            elif len(p_ref_domain) != 0 and len(p_sim_domain) == 0:
                                for i in range(len(p_ref_domain)):
                                    lod.append([p_ref_domain[i][0], p_ref_domain[i][1]])

                            ## In case of both transcripts have more than one domain blocks
                            import pybedtools
                            bed1 = pybedtools.BedTool("\n".join(f"chr1\t{int(start)}\t{int(end)}" for start, end in p_ref_domain), from_string=True)
                            bed2 = pybedtools.BedTool("\n".join(f"chr1\t{int(start)}\t{int(end)}" for start, end in p_sim_domain), from_string=True)
                            
                            lod = []
                            god = []
                            if len(p_ref_domain) >= len(p_sim_domain):
                                for i in range(len(p_ref_domain)):  # Compare length of each block
                                    r_dom = (p_ref_domain[i][1] - p_ref_domain[i][0])
                                    try:
                                        s_dom = (p_sim_domain[i][1] - p_sim_domain[i][0])
                                    except IndexError:
                                        s_dom = 0
                                        
                                    if r_dom > s_dom:
                                        lod.append(i)
                                    elif r_dom < s_dom:
                                        god.append(i)
                                lod_domain = [p_ref_domain[i] for i in lod]
                                god_domain = [p_sim_domain[i] for i in god]
                            
                            else:
                                for i in range(len(p_sim_domain)):  # Compare length of each block
                                    s_dom = (p_sim_domain[i][1] - p_sim_domain[i][0])
                                    try:
                                        r_dom = (p_ref_domain[i][1] - p_ref_domain[i][0])
                                    except IndexError:
                                        r_dom = 0
                                    if r_dom > s_dom:
                                        lod.append(i)
                                    elif r_dom < s_dom:
                                        god.append(i)
                                lod_domain = [p_ref_domain[i] for i in lod]
                                god_domain = [p_sim_domain[i] for i in god]
                            
                            Comp_plot(axs2, lod_domain, "#E72B2B", [n_domain[kinds_domain]])    # Think about dict color map with domain name
                            Comp_plot(axs2, god_domain, "#69EB89", [n_domain[kinds_domain]])    # Think about dict color map with domain name

                            ### Debug
                            # axs3 = plt.subplot(2+(len(n_domain)*2), 1, 3+diff_dom, sharex=ref_axs)
                            # Comp_plot(axs3, p_ref_domain, "#2BE7E7", [n_domain[kinds_domain]])    # Think about dict color map with domain name
                            # axs3 = plt.subplot(2+(len(n_domain)*2), 1, 3+diff_dom+1, sharex=ref_axs)
                            # Comp_plot(axs3, p_sim_domain, "#838383", [n_domain[kinds_domain]])    # Think about dict color map with domain name

                        ## It should be updated, FILTER OUT NO CHANGES AND USING COLOR FOR LoF and GoF
                        if np.abs(np.sum(l_ref_domain) - np.sum(l_sim_domain)) != 0 and \
                            (ref_nmd < 55 and sim_nmd < 55): # It will depict only Domain_alt case
                            diff_dom += 1
                            Domain_plot(diff_dom)

                        fig.suptitle(f'{gene}; Delta L:{round(deltaL,3)}', y=1.1)

                    else:
                        print("refTX does not have domain")
                    
                plt.tight_layout(h_pad=0)
                if not os.path.exists(OUT+"consequence/"):
                    os.mkdir(OUT+"consequence/")
                plt.savefig(OUT+"consequence/{}_comp.pdf".format(gene),
                            bbox_inches="tight")
                print(f"Figure will be saved '{OUT}consequence/'")

        if len(n_domain) > 0:
            Draw_splicemap(sub,sub_ref,sub_sim_tx, "on")

        else:
            print("No mapped functional regions")

Run(args.cano_tx)
# %%
