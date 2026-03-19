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
    parser.add_argument('--splicing_event', '-s', 
                        help='Splicing categories e.g., CA, RI, A3SS, A5SS, MXE', 
                        required=True,
                        type=str)
    parser.add_argument('--gene', '-g', 
                        help='Gene symbol', 
                        required=True,
                        type=str)
    parser.add_argument('--sim_splicing_event', '-sim', 
                        help='Simulated splicing e.g., EI, ES, RI, SI, Ori_A3SS, Alt_A3SS MXE1, MXE2', 
                        required=True,
                        type=str)
    parser.add_argument('--transcript', '-t', 
                        help='Reference transcript that contains same exon structure', 
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

def Run(key, input_gene_name, *sub_splicing):
    """ Draw splicing consequence

    Args:
        key (_type_): Splicing categories (CA, RI, A3SS ...)
        input_gene_name (_type_): Gene_symbol
        sub_splicing (_type_): ES, EI, RI, SI ...

    Returns:
        _type_: _description_
    """
    
    def Load_data(key, remove_info):
        """ Load information data
            It could be smaller (perhaps using cmd?)
        Args:
            key (_type_): Splicing category

        Returns:
            _type_: _description_
        """
        sim_tx = pd.read_csv(args.input+"post_input/merged_w_"+key+".bed",
                                sep="\t", header=None)
        ref = pd.read_csv(args.input+"post_input/merged_wo_"+key+".bed",
                                sep="\t", header=None)
        query = pd.read_csv(args.input+"result/all_{}_Main_table.tsv".format(key),
                            sep="\t", skiprows=1)
        
        query = query[query["5'UTR_difference (nt)"]!="Frame loss"] # Filter out frame loss, because they donot have domain
        query = query[query["ORF"]=="pORF1"]
        query["ORF"] = query["ORF"].apply(lambda x : 3 if x =="pORF1" else (2 if x == "pORF2" else 1))
        query = query[["LongID","Gene symbol","Simulated_event","AUG (Ref-Sim)","Stop (Ref-Sim)","CDS_difference (nt)","Reference_transcript","ORF","Domain_change_rate","Functional_class","Probability_of_NMD"]]
        query["simStop"] = query["Stop (Ref-Sim)"].str.split("-").str[1]
        query["Stop (Ref-Sim)"] = query["Stop (Ref-Sim)"].str.split("-").str[0]
        query["simStart"] = query["AUG (Ref-Sim)"].str.split("-").str[1]
        query["AUG (Ref-Sim)"] = query["AUG (Ref-Sim)"].str.split("-").str[0]
        query.columns = ["ID","Gene symbol","event_type","Start","STOP","dAA","ENST","ORF","delta L","DOA_types","pNMD","simSTOP","simStart"]
        query["SID"] = query["ID"]
        query["ID"] = query["ID"]+"|"+query["ENST"]+"|"+query["Start"]+"|"+query["simStart"]
        query[["simSTOP","STOP","Start","simStart"]] = query[["simSTOP","STOP","Start","simStart"]].astype(float)

        w_pfam = pd.read_csv(args.input+"table/{}_w_Pfam.txt".format(key),
                        sep="\t", header=None)
        wo_pfam = pd.read_csv(args.input+"table/{}_wo_Pfam.txt".format(key),
                        sep="\t", header=None)
        
        w_pfam = w_pfam[~w_pfam[4].isin(remove_info)]
        wo_pfam = wo_pfam[~wo_pfam[4].isin(remove_info)]

        w_pfam[3] = w_pfam[3] +" ("+w_pfam[4]+")"
        wo_pfam[3] = wo_pfam[3] +" ("+wo_pfam[4]+")"

        return sim_tx, ref, query, w_pfam, wo_pfam


    sim_tx, ref, query, w_pfam, wo_pfam = Load_data(key, args.remove_info)

    ### Specify your input
    if len(sub_splicing) > 0 and len(sub_splicing) < 2:
        query = query[query["event_type"]==sub_splicing[0]]
    elif len(sub_splicing) > 1:
        query = query[query["ENST"]==sub_splicing[1]]
        query = query[query["event_type"]==sub_splicing[0]]
    else:
        pass
    query["gene"] = query["ID"].str.split(";").str[1]

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
                ########################################################################
                if sub_splicing[0] in ["Alt_A3SS", "Alt_A5SS", "RI"]:
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#FFFFFF",
                                            zorder=5, linestyle=":", linewidth=1.5,
                                            hatch="///")
                        
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="none",
                                            zorder=4, linestyle="-", linewidth=1.5)
                
                elif sub_splicing[0] in ["Ori_A3SS", "Ori_A5SS", "SI"]:
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='none',
                                            zorder=5, linestyle="-", linewidth=1.5)
                    
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor="#FFFFFF",
                                            zorder=5, linestyle=":", linewidth=1.5,
                                            hatch="///")
                
                elif sub_splicing[0].startswith("MXE"):
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        if domain[-1] == "on":
                            fc = "#FFFFFF"
                            canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                                edgecolor="#000000", facecolor=fc,
                                                zorder=5, linestyle=":", linewidth=1.5,
                                                hatch="///")
                        else:
                            fc = "none"
                            canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                                edgecolor="#000000", facecolor=fc,
                                                zorder=4, linestyle="-", linewidth=1.5)
                    
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        if domain[-1] == "on":
                            fc = "none"
                            canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                                edgecolor="#000000", facecolor=fc,
                                                zorder=4, linestyle="-", linewidth=1.5)
                        else:
                            fc = "#FFFFFF"
                            canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                                edgecolor="#000000", facecolor=fc,
                                                zorder=5, linestyle=":", linewidth=1.5,
                                                hatch="///")

                elif sub_splicing[0] == "ES":
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='none',
                                            zorder=4, linestyle="-", linewidth=1.5)
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='#FFFFFF',
                                            zorder=5, linestyle=":", linewidth=1.5,
                                            hatch="///")
                
                elif sub_splicing[0] == "EI":
                    if len(domain) > 2 and exon_color != "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='#FFFFFF',
                                            zorder=5, linestyle=":", linewidth=1.5,
                                            hatch="///")
                    elif len(domain) > 2 and exon_color == "#73ABFF":
                        canvas.fill_between(span, ylims['exon_min'], ylims['exon_max'],
                                            edgecolor="#000000", facecolor='none',
                                            zorder=4, linestyle="-", linewidth=1.5)
                ########################################################################
                
                ########################################################################
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
                    
                    if span[0] == 0:    # Consequence annotation
                        value_dict = {"PTC removal":"NMD",
                                      "NMD":"Intact",
                                      "LoD":"GoD",
                                      "GoD":"LoD",
                                      "CDS_alt":"CDS_alt",
                                      "UTR_alt":"UTR_alt",
                                      "no_changes":"No Changes"} # BUG fix 25.01.27
                        if exon_color != "#73ABFF": # Ref-TX
                            canvas.text(.96,.2,f"{value_dict[doa_types]}",
                                        fontsize=10, zorder=5,
                                        transform=canvas.transAxes)
                        else:
                            canvas.text(.96,.2,f"{doa_types}",
                                    fontsize=10, zorder=5,
                                    transform=canvas.transAxes)
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
                mxe_exon = 0    # For MXE exon marking
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

                    ## Drawing the exon structure with changes
                    elif len(pos_info) >= 3:    # MXE, EI, A3/5SS
                        if sub_splicing[0] == ("MXE1"):   # MXE cases
                            if i in [pos_info[2], pos_info[2]+1]:    # Changed position
                                new_pos_info = [j for j in pos_info]
                                if mxe_exon == 0:   # For 1st MXE
                                    new_pos_info.append("on")
                                else:   # For 2nd MXE
                                    new_pos_info.append("off")
                                draw_exon(exonpos[i], canvas, *new_pos_info)
                                mxe_exon += 1
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1]) # It only has two pos_info (==len(domain)==2 in draw_exon function)
                        
                        elif sub_splicing[0] == ("MXE2"):   # MXE cases
                            if i in [pos_info[2], pos_info[2]+1]:    # Changed position
                                new_pos_info = [j for j in pos_info]
                                if mxe_exon == 0:   # For 1st MXE
                                    new_pos_info.append("off")
                                else:   # For 2nd MXE
                                    new_pos_info.append("on")
                                draw_exon(exonpos[i], canvas, *new_pos_info)
                                mxe_exon += 1
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1]) # It only has two pos_info (==len(domain)==2 in draw_exon function)
                        
                        elif sub_splicing[0] in ["EI","ES"]:   # For EI/S cases
                            if i == pos_info[2]:  # It has changed exon
                                draw_exon(exonpos[i], canvas, *pos_info)
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1]) # It only has two pos_info (==len(domain)==2 in draw_exon function)
                        
                        elif sub_splicing[0].endswith("A3SS"):  # For A3/5SS or RI/SI cases
                            if i == pos_info[2]:  # It has changed exon
                                try:    # Draw 
                                    draw_exon(pos_info[3], canvas, *pos_info)   # AA splice block
                                    draw_exon(pos_info[4], canvas, pos_info[0], pos_info[1])   # after AA splice block
                                except:
                                    pass
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1]) # It only has two pos_info (==len(domain)==2 in draw_exon function)
                        
                        elif sub_splicing[0].endswith("A5SS") or \
                            sub_splicing[0].endswith("I"):   # For A5SS or RI/SI cases
                            if i == pos_info[2]:  # It has changed exon
                                try:    # Draw 
                                    draw_exon(pos_info[3], canvas, pos_info[0], pos_info[1])   # before splice block, If "pos_info[0] and [1] are replaced by *pos_info, it have dotted outline"
                                    draw_exon(pos_info[4], canvas, *pos_info)   # AA splice block
                                    draw_exon(pos_info[5], canvas, pos_info[0], pos_info[1])   # after AA splice block
                                except:
                                    pass
                            else:   # If it is not changed exon
                                draw_exon(exonpos[i], canvas, pos_info[0], pos_info[1]) # It only has two pos_info (==len(domain)==2 in draw_exon function)

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
        ##############################
        strand = final[5].values[0]
        for line in range(final.shape[0]):  ## Individual transcript
            if line > 0:    ## From 2nd exon
                if strand == "+":
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                  (pfam[9]==final.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                    ## Previous exon end - current exon start
                    intron = np.abs(int(final.iloc[line-1,2] - int(final.iloc[line,1])))
                    diff.append(intron)
                    ## Subtract cumulative intron from each exon coordinate
                    ## zero is the coordinate of start of first exon
                    new_line.append([int(final.iloc[line,1])-np.sum(diff)-zero, 
                                     int(final.iloc[line,2])-np.sum(diff)-zero])
                    ## Debug
                    # print(final.iloc[line,1],final.iloc[line,2])
                    # print(int(final.iloc[line,1])-np.sum(diff)-zero,int(final.iloc[line,2])-np.sum(diff)-zero)
                    ###########

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
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                  (pfam[9]==final.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                    intron = np.abs(np.abs(int(final.iloc[line-1,1])-zero) - np.abs(int(final.iloc[line,2])-zero))
                    diff.append(intron)
                    new_line.append([np.abs(int(final.iloc[line,2])-zero)-np.sum(diff), 
                                     np.abs(int(final.iloc[line,1])-zero)-np.sum(diff)])
                    
                    ## Debug
                    # print(final.iloc[line,2],final.iloc[line,1])
                    # print(np.abs(int(final.iloc[line,2])-zero)-np.sum(diff),np.abs(int(final.iloc[line,1])-zero)-np.sum(diff))
                    ###########
                    
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

            else:   ## For 1st exon
                pre_domain = " "
                if strand == "+":
                    zero = int(final.iloc[line,1])  # The most left exon's start site for align in figure
                    new_line.append([int(final.iloc[line,1])-zero,
                                    int(final.iloc[line,2])-zero])
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                  (pfam[9]==final.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                    
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
                    zero = int(final.iloc[line,2])  # The most right exon's end site
                    new_line.append([np.abs(int(final.iloc[line,2])-zero),
                                    np.abs(int(final.iloc[line,1])-zero)])
                    domain = pfam[(pfam[8]==final.iloc[line,1]) &
                                  (pfam[9]==final.iloc[line,2])][[0,1,2,3,8,9]].drop_duplicates()
                    
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


    test_list = query.iloc[:,0].unique()
    # sid = test_list[0].split("|")[0]
    test_list = query[query["Gene symbol"] == input_gene_name]["ID"].unique()

    if len(test_list) == 0:
        print("\nERROR: Please make sure your Gene symbol tab in the Main_table.tsv!!\n")
        sys.exit(1)
    
    ################################ TEMP
    ## Remove Biding site in figures, it spends lots of space
    ## Capture target domains which are contained in target transcript
    w_pfam = w_pfam[(w_pfam[14].isin(query["SID"].unique())) &
                    (w_pfam[11]==sub_splicing[1]) &
                    (w_pfam[13]==sub_splicing[0]) &
                    (w_pfam[4]!="binding")]
    wo_pfam = wo_pfam[(wo_pfam[14].isin(query["SID"].unique())) &
                      (wo_pfam[11]==sub_splicing[1]) &
                      (wo_pfam[13]==sub_splicing[0]) &
                      (wo_pfam[4]!="binding")]
    n_domain = wo_pfam[3].unique()
    ################################

    OUT = args.input+"figure/"
    for id in test_list:    # Certain splicing event with one of the possible ORF
        LID = id.split("|")[0]
        sub = query[query["ID"]==id]

        ######## Functional Annotation
        doa_types = str(sub["DOA_types"].values[0])
        if sub["pNMD"].values[0] == "-1":
            doa_types = "NMD"
        elif sub["pNMD"].values[0] == "1":
            doa_types = "PTC removal"
        ########

        deltaL = float(sub["delta L"].values[0])
        sub_sim_tx = sim_tx[(sim_tx[4]==sub["ENST"].values[0]) &
                            (sim_tx[6]==sub["event_type"].values[0]) &
                            (sim_tx[7]==sub["SID"].values[0])].copy() # ASS events can be distinguished by event type, but other event need ID
        sub_ref = ref[(ref[4]==sub["ENST"].values[0]) & 
                      (ref[6]==sub["event_type"].values[0])&
                      (ref[7]==sub["SID"].values[0])].copy()

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
                    gene = fig_input["Gene symbol"].values[0]
                    tx = fig_input["ENST"].values[0]
                    eventid = fig_input["event_type"].values[0]

                    ## Make an input (bed) for visualization
                    ref_bed, ref_domain, ref_dname = Make_query(ref, wo_pfam[(wo_pfam[3]==n_domain[kinds_domain]) & (wo_pfam[14].isin(sub["SID"].unique()))])
                    sim_bed, sim_domain, sim_dname = Make_query(sim_tx, w_pfam[(w_pfam[3]==n_domain[kinds_domain]) & (w_pfam[14].isin(sub["SID"].unique()))])
                    # print(n_domain[kinds_domain])
                    # print(ref_bed)
                    # print(ref_domain)
                    # print(sim_bed)
                    # print(sim_domain)

                    ## Check different exon block
                    if len(ref_bed) > 0:    # Markup changed regions
                        Make_fig = 0
                        for r,s in zip(ref_bed,sim_bed):
                            if r[0] != s[0] or \
                                r[1] != s[1]:
                                event_pos_ref = ref_bed.index(r)    # Position of DS on Ref TX
                                event_pos_query = sim_bed.index(s)  # Position of DS on Sim TX
                                Make_fig += 1   # Make altered regions due to the given DS event
                                break   # there is at least one different CDS alterations)
                        
                        ################################
                        ## Make a prunned bed
                        ## Make a dataframe for subsequent analysis
                        def Prunning(domain_bed, dname, st, end, diff_exon, space, event):
                            ## Add truncated domain
                            ## Motif could be located on UTR regions
                            l_domain = []   # length
                            p_domain_bed = []   # domain position
                            p_dname = []    # domain name
                            for index in range(len(domain_bed)):
                                if not event.startswith('MXE'):
                                    ## Adjust position of Protein in figure
                                    if int(domain_bed[index][0]) >= int(diff_exon[0]):   # if domain is located after diff exon
                                        domain_bed[index][0] += space
                                        domain_bed[index][1] += space
                                    elif int(domain_bed[index][0]) < int(diff_exon[0]) and \
                                        int(domain_bed[index][1]) > int(diff_exon[0]):    # partial overlap
                                        domain_bed[index][0] = int(diff_exon[0]) + int(space)   # domain_bed[index][0]~diff_exon[0] were not included
                                        domain_bed[index][1] += space

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

                    ###################################
                    ## Make an updated bed to draw skipped/included parts
                    ########
                        # print(event_pos_ref, event_pos_query)
                        space = np.abs(float(LID.split(";")[6])-float(LID.split(";")[7]))+1
                        if fig_input["event_type"].unique()[0] in ["EI","Alt_A3SS"]:

                            if fig_input["Start"].values[0] >= sim_bed[event_pos_query][0]:
                                adj_ref_start = fig_input["Start"].values[0]+space
                                adj_ref_stop = fig_input["STOP"].values[0]+space
                            elif fig_input["Start"].values[0] < sim_bed[event_pos_query][0] and\
                                 fig_input["STOP"].values[0] >= sim_bed[event_pos_query][0]:
                                adj_ref_start = fig_input["Start"].values[0]
                                adj_ref_stop = fig_input["STOP"].values[0]+space
                            else:
                                adj_ref_start = fig_input["Start"].values[0]
                                adj_ref_stop = fig_input["STOP"].values[0]
                            adj_sim_start = fig_input["simStart"].values[0]
                            adj_sim_stop = fig_input["simSTOP"].values[0]

                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               adj_ref_start, 
                                                                               adj_ref_stop,
                                                                               ref_bed[event_pos_ref],
                                                                               space,
                                                                               fig_input["event_type"].unique()[0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               adj_sim_start, 
                                                                               adj_sim_stop,
                                                                               [0,0],
                                                                               0,
                                                                               fig_input["event_type"].unique()[0])
                            
                        elif fig_input["event_type"].unique()[0] in ["ES","Ori_A3SS"]:

                            if fig_input["simStart"].values[0] >= ref_bed[event_pos_ref][0]:
                                adj_sim_start = fig_input["simStart"].values[0]+space
                                adj_sim_stop = fig_input["simSTOP"].values[0]+space
                            elif fig_input["simStart"].values[0] < ref_bed[event_pos_ref][0] and\
                                 fig_input["simSTOP"].values[0] >= ref_bed[event_pos_ref][0]:
                                adj_sim_start = fig_input["simStart"].values[0]
                                adj_sim_stop = fig_input["simSTOP"].values[0]+space
                            else:
                                adj_sim_start = fig_input["simStart"].values[0]
                                adj_sim_stop = fig_input["simSTOP"].values[0]
                            adj_ref_start = fig_input["Start"].values[0]
                            adj_ref_stop = fig_input["STOP"].values[0]

                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               adj_ref_start, 
                                                                               adj_ref_stop,
                                                                               [0,0],
                                                                               0,
                                                                               fig_input["event_type"].unique()[0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               adj_sim_start, 
                                                                               adj_sim_stop,
                                                                               sim_bed[event_pos_query],
                                                                               space,
                                                                               fig_input["event_type"].unique()[0])
                            
                        elif fig_input["event_type"].unique()[0] in ["Alt_A5SS","RI"]:

                            if fig_input["Start"].values[0] >= ref_bed[event_pos_ref][1]:
                                adj_ref_start = fig_input["Start"].values[0]+space
                                adj_ref_stop = fig_input["STOP"].values[0]+space
                            elif fig_input["Start"].values[0] < ref_bed[event_pos_ref][1] and\
                                 fig_input["STOP"].values[0] >= ref_bed[event_pos_ref][1]:
                                adj_ref_start = fig_input["Start"].values[0]
                                adj_ref_stop = fig_input["STOP"].values[0]+space
                            else:
                                adj_ref_start = fig_input["Start"].values[0]
                                adj_ref_stop = fig_input["STOP"].values[0]
                            adj_sim_start = fig_input["simStart"].values[0]
                            adj_sim_stop = fig_input["simSTOP"].values[0]

                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               adj_ref_start, 
                                                                               adj_ref_stop,
                                                                               [ref_bed[event_pos_ref][1]],
                                                                               space,
                                                                               fig_input["event_type"].unique()[0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               adj_sim_start, 
                                                                               adj_sim_stop,
                                                                               [0,0],
                                                                               0,
                                                                               fig_input["event_type"].unique()[0])
 
                        elif fig_input["event_type"].unique()[0] in ["Ori_A5SS","SI"]:

                            if fig_input["simStart"].values[0] >= ref_bed[event_pos_ref][0]:
                                adj_sim_start = fig_input["simStart"].values[0]+space
                                adj_sim_stop = fig_input["simSTOP"].values[0]+space
                            elif fig_input["simStart"].values[0] < ref_bed[event_pos_ref][0] and\
                                 fig_input["simSTOP"].values[0] >= ref_bed[event_pos_ref][0]:
                                adj_sim_start = fig_input["simStart"].values[0]
                                adj_sim_stop = fig_input["simSTOP"].values[0]+space
                            else:
                                adj_sim_start = fig_input["simStart"].values[0]
                                adj_sim_stop = fig_input["simSTOP"].values[0]
                            adj_ref_start = fig_input["Start"].values[0]
                            adj_ref_stop = fig_input["STOP"].values[0]

                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(ref_domain, ref_dname,
                                                                               adj_ref_start, 
                                                                               adj_ref_stop,
                                                                               [0,0],
                                                                               0,
                                                                               fig_input["event_type"].unique()[0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(sim_domain, sim_dname,
                                                                               adj_sim_start, 
                                                                               adj_sim_stop,
                                                                               [sim_bed[event_pos_query][1]],
                                                                               space,
                                                                               fig_input["event_type"].unique()[0])

                        else:   # MXE 1 or 2
                            if fig_input["event_type"].unique()[0] == "MXE2":
                                mxe1_leng = np.abs(ref_bed[event_pos_ref][0] - ref_bed[event_pos_ref][1])
                                mxe2_leng = np.abs(sim_bed[event_pos_query][0] - sim_bed[event_pos_query][1])
                                space1 = mxe1_leng  # ref block
                                space2 = mxe2_leng  # sim block
                                end_of_MXE1 = ref_bed[event_pos_ref][1]
                                mxe_bed = ref_bed[:event_pos_ref+1] + \
                                        [[sim_bed[event_pos_query][0]+space1, sim_bed[event_pos_query][1]+space1]] + \
                                        [(i+space2) for i in ref_bed[event_pos_ref+1:]] # Add all MXE to draw same length of Ref- and Sim-TX

                                if fig_input["simStart"].values[0] >= ref_bed[event_pos_ref][0]:
                                    adj_sim_start = fig_input["simStart"].values[0]+space1
                                    adj_sim_stop = fig_input["simSTOP"].values[0]+space1
                                elif fig_input["simStart"].values[0] < ref_bed[event_pos_ref][0] and\
                                    fig_input["simSTOP"].values[0] >= ref_bed[event_pos_ref][0]:
                                    adj_sim_start = fig_input["simStart"].values[0]
                                    adj_sim_stop = fig_input["simSTOP"].values[0]+space1
                                else:
                                    adj_sim_start = fig_input["simStart"].values[0]
                                    adj_sim_stop = fig_input["simSTOP"].values[0]
                                
                                if fig_input["Start"].values[0] >= sim_bed[event_pos_query][0]:
                                    adj_ref_start = fig_input["Start"].values[0]+space2
                                    adj_ref_stop = fig_input["STOP"].values[0]+space2
                                elif fig_input["Start"].values[0] < sim_bed[event_pos_query][0] and\
                                    fig_input["STOP"].values[0] >= sim_bed[event_pos_query][0]:
                                    adj_ref_start = fig_input["Start"].values[0]
                                    adj_ref_stop = fig_input["STOP"].values[0]+space2
                                else:
                                    adj_ref_start = fig_input["Start"].values[0]
                                    adj_ref_stop = fig_input["STOP"].values[0]

                            if fig_input["event_type"].unique()[0] == "MXE1":
                                mxe1_leng = np.abs(sim_bed[event_pos_query][0] - sim_bed[event_pos_query][1])
                                mxe2_leng = np.abs(ref_bed[event_pos_ref][0] - ref_bed[event_pos_ref][1])
                                space1 = mxe2_leng  # ref block
                                space2 = mxe1_leng  # sim block
                                end_of_MXE1 = sim_bed[event_pos_query][1]
                                mxe_bed = ref_bed[:event_pos_ref] + \
                                        [[sim_bed[event_pos_query][0],sim_bed[event_pos_query][1]]] + \
                                        [(i+space2) for i in ref_bed[event_pos_ref:]]
                                
                                if fig_input["simStart"].values[0] >= ref_bed[event_pos_ref][0]:
                                    adj_sim_start = fig_input["simStart"].values[0]+space1
                                    adj_sim_stop = fig_input["simSTOP"].values[0]+space1
                                elif fig_input["simStart"].values[0] < ref_bed[event_pos_ref][0] and\
                                    fig_input["simSTOP"].values[0] >= ref_bed[event_pos_ref][1]:
                                    adj_sim_start = fig_input["simStart"].values[0]
                                    adj_sim_stop = fig_input["simSTOP"].values[0]+space1
                                else:
                                    adj_sim_start = fig_input["simStart"].values[0]
                                    adj_sim_stop = fig_input["simSTOP"].values[0]
                                
                                if fig_input["Start"].values[0] >= sim_bed[event_pos_query][0]:
                                    adj_ref_start = fig_input["Start"].values[0]+space2
                                    adj_ref_stop = fig_input["STOP"].values[0]+space2
                                elif fig_input["Start"].values[0] < sim_bed[event_pos_query][0] and\
                                    fig_input["STOP"].values[0] >= sim_bed[event_pos_query][1]:
                                    adj_ref_start = fig_input["Start"].values[0]
                                    adj_ref_stop = fig_input["STOP"].values[0]+space2
                                else:
                                    adj_ref_start = fig_input["Start"].values[0]
                                    adj_ref_stop = fig_input["STOP"].values[0]
                            
                            ## Set the domain position to compatible with the mxe_bed
                            updated_sim_domain = []
                            for uncorr in sim_domain:
                                if uncorr[0] >= end_of_MXE1:
                                    corr = list(uncorr+mxe2_leng)
                                    updated_sim_domain.append(corr)
                                elif uncorr[0] < end_of_MXE1 and uncorr[1] > end_of_MXE1:
                                    corr1 = [uncorr[0], end_of_MXE1]
                                    corr2 = [end_of_MXE1 + mxe2_leng, uncorr[1] + mxe2_leng]
                                    updated_sim_domain.append(corr1)
                                    updated_sim_domain.append(corr2)
                                else:
                                    updated_sim_domain.append(uncorr)
                            
                            updated_ref_domain = []
                            for uncorr in ref_domain:
                                if uncorr[0] >= end_of_MXE1:
                                    corr = list(uncorr+mxe1_leng)
                                    updated_ref_domain.append(corr)
                                elif uncorr[0] <= end_of_MXE1 and uncorr[1] > end_of_MXE1:
                                    corr1 = [uncorr[0], end_of_MXE1]
                                    corr2 = [end_of_MXE1 + mxe1_leng, uncorr[1] + mxe1_leng]
                                    updated_ref_domain.append(corr1)
                                    updated_ref_domain.append(corr2)
                                else:
                                    updated_ref_domain.append(uncorr)
                                    
                            p_ref_domain, p_ref_dname, l_ref_domain = Prunning(updated_ref_domain, ref_dname,
                                                                               adj_ref_start, 
                                                                               adj_ref_stop,
                                                                               sim_bed[event_pos_query],
                                                                               space2,
                                                                               fig_input["event_type"].unique()[0])
                            p_sim_domain, p_sim_dname, l_sim_domain = Prunning(updated_sim_domain, sim_dname,
                                                                               adj_sim_start, 
                                                                               adj_sim_stop,
                                                                               ref_bed[event_pos_ref],
                                                                               space1,
                                                                               fig_input["event_type"].unique()[0])
                          
                    ###################################
                        if (len(p_ref_domain) > 0 or len(p_sim_domain) > 0) and fig_on == "on":
                            if n_row == 0:
                                ## Main complex figure
                                ref_axs = plt.subplot(2+(len(n_domain)*2), 1, 1)
                                if Make_fig == 1:
                                    ## If sim splicing is inclusion form, ref tx should have same tx structure
                                    ## Draw Ref-TX
                                    if fig_input["event_type"].unique()[0] in ["EI"]:
                                        ## Comp_plot: line 123
                                        Comp_plot(ref_axs, sim_bed, "#FF7878", 
                                                  adj_ref_start, 
                                                  adj_ref_stop, 
                                                  event_pos_ref)

                                    elif fig_input["event_type"].unique()[0] in ["ES"]:
                                        Comp_plot(ref_axs, ref_bed, "#FF7878", 
                                                  adj_ref_start, 
                                                  adj_ref_stop, 
                                                  event_pos_ref)
                                        
                                    elif fig_input["event_type"].unique()[0] in ["Alt_A5SS","Alt_A3SS"]:   # It should be tested
                                        Comp_plot(ref_axs, sim_bed, "#FF7878", 
                                                  adj_ref_start, 
                                                  adj_ref_stop, 
                                                  event_pos_ref,
                                                  [sim_bed[event_pos_query][0], sim_bed[event_pos_query][0]+space],
                                                  [sim_bed[event_pos_query][0]+space, sim_bed[event_pos_query][1]])

                                    elif fig_input["event_type"].unique()[0] in ["Ori_A5SS","Ori_A3SS"]:
                                        Comp_plot(ref_axs, ref_bed, "#FF7878", 
                                                  adj_ref_start, 
                                                  adj_ref_stop, 
                                                  event_pos_ref,
                                                  [sim_bed[event_pos_query][0], sim_bed[event_pos_query][0]+space],
                                                  [sim_bed[event_pos_query][0]+space, ref_bed[event_pos_ref][1]])

                                    elif fig_input["event_type"].unique()[0] in ["RI","SI"]:   # RI and SI
                                        if fig_input["event_type"].unique()[0] == "RI":
                                            Comp_plot(ref_axs, sim_bed, "#FF7878", 
                                                      adj_ref_start, 
                                                      adj_ref_stop, 
                                                      event_pos_ref,
                                                      [ref_bed[event_pos_ref][0], ref_bed[event_pos_ref][1]],
                                                      [ref_bed[event_pos_ref][1], ref_bed[event_pos_ref][1]+space],
                                                      [ref_bed[event_pos_ref][1]+space, sim_bed[event_pos_query][1]])
                                            
                                        else:
                                            Comp_plot(ref_axs, ref_bed, "#FF7878", 
                                                      adj_ref_start, 
                                                      adj_ref_stop, 
                                                      event_pos_ref,
                                                      [sim_bed[event_pos_query][0], sim_bed[event_pos_query][1]],
                                                      [sim_bed[event_pos_query][1], sim_bed[event_pos_query][1]+space],
                                                      [sim_bed[event_pos_query][1]+space, ref_bed[event_pos_ref][1]])

                                    else:   # MXE
                                        Comp_plot(ref_axs, mxe_bed, "#FF7878", 
                                                  adj_ref_start, 
                                                  adj_ref_stop, 
                                                  event_pos_ref)

                                plt.text(-0.03,0.2,"{}".format("Ref-Tx"), 
                                        transform=ref_axs.transAxes,
                                        fontsize=10)
                                
                                ## Draw Sim-TX
                                sim_axs = plt.subplot(2+(len(n_domain)*2), 1, 2, sharex=ref_axs)
                                if fig_input["event_type"].unique()[0] in ["ES"]:
                                    Comp_plot(sim_axs, ref_bed, "#73ABFF", 
                                              adj_sim_start, 
                                              adj_sim_stop, 
                                              event_pos_query)
                                elif fig_input["event_type"].unique()[0] in ["Ori_A3SS", "Ori_A5SS"]:
                                    Comp_plot(sim_axs, ref_bed, "#73ABFF", 
                                              adj_sim_start, 
                                              adj_sim_stop, 
                                              event_pos_ref,
                                              [sim_bed[event_pos_query][0], sim_bed[event_pos_query][0]+space],
                                              [sim_bed[event_pos_query][0]+space, ref_bed[event_pos_ref][1]])
                                    
                                elif fig_input["event_type"].unique()[0] in ["EI"]:
                                    Comp_plot(sim_axs, sim_bed, "#73ABFF", 
                                              adj_sim_start, 
                                              adj_sim_stop, 
                                              event_pos_query)
                                
                                elif fig_input["event_type"].unique()[0] in ["RI","SI"]:
                                    if fig_input["event_type"].unique()[0] == "RI":
                                        Comp_plot(sim_axs, sim_bed, "#73ABFF", 
                                                  adj_sim_start, 
                                                  adj_sim_stop, 
                                                  event_pos_query,
                                                  [ref_bed[event_pos_ref][0], ref_bed[event_pos_ref][1]],
                                                  [ref_bed[event_pos_ref][1], ref_bed[event_pos_ref][1]+space],
                                                  [ref_bed[event_pos_ref][1]+space, sim_bed[event_pos_query][1]])
                                                
                                    else:
                                        Comp_plot(sim_axs, ref_bed, "#73ABFF", 
                                                  adj_sim_start, 
                                                  adj_sim_stop, 
                                                  event_pos_query,
                                                  [sim_bed[event_pos_query][0], sim_bed[event_pos_query][1]],
                                                  [sim_bed[event_pos_query][1], sim_bed[event_pos_query][1]+space],
                                                  [sim_bed[event_pos_query][1]+space, ref_bed[event_pos_ref][1]])
                                
                                elif fig_input["event_type"].unique()[0] in ["Alt_A3SS","Alt_A5SS"]:
                                    Comp_plot(sim_axs, sim_bed, "#73ABFF", 
                                              adj_sim_start, 
                                              adj_sim_stop, 
                                              event_pos_query,
                                              [sim_bed[event_pos_query][0], sim_bed[event_pos_query][0]+space],
                                              [sim_bed[event_pos_query][0]+space, sim_bed[event_pos_query][1]])
                                
                                else:   # MXE
                                    Comp_plot(sim_axs, mxe_bed, "#73ABFF", 
                                              adj_sim_start, 
                                              adj_sim_stop, 
                                              event_pos_query)
                                    
                                plt.text(-0.03,0.2,"{}".format("Sim-Tx"), 
                                        transform=sim_axs.transAxes,
                                        fontsize=10)
                                n_row += 1  # For next domain, tx drawing function will be skipped


                            ## Draw protein domain
                            ## Last exon-exon junction - stop codon, it may be less than 55 it will be considered as NMD
                            ref_nmd = ref_bed[-1][0] - fig_input["STOP"].values[0]
                            sim_nmd = sim_bed[-1][0] - fig_input["simSTOP"].values[0]
                            
                            def Domain_plot(diff_dom):
                                axs2 = plt.subplot(2+(len(n_domain)*2), 1, 3+diff_dom, sharex=ref_axs)
                                ## Check here, domain block concordance
                                ## It should be considered skipped or included exon's length. It is key!!!!
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

                                # shared
                                # overlap = bed1.intersect(bed2)

                                # ref-sim
                                loss_of_domain = bed1.subtract(bed2)
                                
                                # sim-ref
                                gain_of_domain = bed2.subtract(bed1)
                            
                                for interval in loss_of_domain:
                                    lod.append([float(interval[1]), float(interval[2])])
                                
                                for interval in gain_of_domain:
                                    god.append([float(interval[1]), float(interval[2])])


                                Comp_plot(axs2, lod, "#E72B2B", [n_domain[kinds_domain]])    # Think about dict color map with domain name
                                Comp_plot(axs2, god, "#69EB89", [n_domain[kinds_domain]])    # Think about dict color map with domain name

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

                            fig.suptitle(f'{gene}; {tx}; sim{eventid}; Delta L:{round(deltaL,3)}\n{LID}', y=1.)

                        else:
                            print("refTX does not have domain")
                    else:
                        plt.close() # Temporal
                
                plt.tight_layout(h_pad=0)
                if not os.path.exists(OUT+"consequence/"):
                    os.mkdir(OUT+"consequence/")

                if len(sub_splicing) > 1:
                    plt.savefig(OUT+"consequence/{}_{}_splicing_map.pdf".format(fig_input["SID"].values[0],sub_splicing[1]),
                                bbox_inches="tight")
                    print(f"Figure will be saved '{OUT}consequence/'")
                elif len(sub_splicing) == 1:
                    plt.savefig(OUT+"consequence/{}_{}_splicing_map.pdf".format(fig_input["SID"].values[0],sub_splicing[0]),
                                bbox_inches="tight")
                    print(f"Figure will be saved '{OUT}consequence/'")
                else:
                    plt.savefig(OUT+"consequence/{}_splicing_map.pdf".format(fig_input["ID"].values[0]),
                                bbox_inches="tight")
                    print(f"Figure will be saved '{OUT}consequence/'")

        if len(n_domain) > 0:
            Draw_splicemap(sub,sub_ref,sub_sim_tx, "on")
        else:
            print("No mapped functional regions")


Run(args.splicing_event, args.gene, args.sim_splicing_event, args.transcript)
