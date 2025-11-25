## What is the SpliceDecoder?
<img width="2436" height="1359" alt="image" src="https://github.com/user-attachments/assets/6d2ba86c-181c-41da-aa0f-a6205cf87509" />

* Splice decoder provides functional annotation for your differential splicing events (DSEs)
* The functional annotation contains NMD probability, alterations in functional domains (such as DNA binding, motif, regions, protein domain, and so on), CDS/UTR alterations, and effect score
* The effect score can be used to prioritize and choose the most representative functional consequences of your DSEs
* Currently, SpliceDecoder supports **hg38** and **mm10** genome

## Workflow overview
<img width="1407" height="509" alt="image" src="https://github.com/user-attachments/assets/ef6b3535-982c-483c-b30f-d2f486092c4d" />


1. Generate All Possible Splicing Cases (Processing input) : This step makes proper format of input data from the output of event-based splicing tools
2. Map Splicing Cases (Mapping DSEs and ORFs): This step explores the given transcriptome (.GTF) to find Ref-TX (Reference transcript, it contains perfectly matched exon structure for the given DSE) and assign the best three open reading frames (ORFs)
3. Simulate Splicing Events (Simulation): Based on the Ref-TXs and their ORFs, this step perform simulation of alternative splicing (e.g., if the Ref-TX has exon inclusion (EI) form, this step makes a simulated transcript (Sim-TX) with exon skipped (ES) form)
4. Functional Annotation (Annotation): Based on the Uniprot DB, SpliceDecoder assigns known functional domains and estimates functional changes between Ref-TX and Sim-TX
5. DSEs with Effect Score (Scoring): SpliceDecoder assigns an effect score to each DSE based on multiple biological factors, enabling prioritization of your DSEs

<br>

## Quick start (conda is required)
* SpliceDecoder can be downloaded from https://github.com/hyeon9/SpliceDecoder/
* Install SpliceDecoder by using the install script
  
      cd ./SpliceDecoder && bash install.sh

* To perform a test run, you can use the provided toy_data
* You can make a toy configuration file through an interactive way [You can find more details here](#guide-for-making-config-file)

      cd code/
      bash Make_config.sh toy

* If you successfully created `Your_toy.config`, you can run SpliceDecoder
* The steps are intended to be executed in order, so it is recommended to use `all`

      bash Main.sh all ${Your_toy.config}

* If your test run with the toy data finishes successfully, [you will see the following output files](#outputs) (except for the Effect_score.tsv)
* To run with your own data, create a configuration file using the `Make_config.sh` or modifying `example.config`

      cd code/
      bash Make_config.sh
      bash Main.sh all ${Your.config}

      OR

      vi example.config
      mv example.config ${Your.config}

* Then, use this command to submit your job if you are using SLURM

      sbatch Main.sh {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring | all} ${Your.config}


* If needed, you can run a specific step by selecting one of the following: `Make_input`, `DS_mapping`, `ORF_mapping`, `Simulation` and `Scoring`

      bash Main.sh {Make_input | DS_mapping | ORF_mapping | Simulation | Scoring} ${Your.config}

* If you want to annotate a transcript-centric data [you can find more details here](./code/transcript-toolkit/README.md)

<br>

## Guide for making config file
* Make_config.sh will ask..

      ? Specify your config file name (e.g. HGjob)
      > You just need to specify your config file

      ? Enter the path of SpliceDecoder (e.g. /User/usr/Tool/SpliceDecoder-main/)
      > You just need to specify the install path of SpliceDecoder
  
      ? Enter your working directory (e.g. /User/usr/Tool/SpliceDecoder-main/project1)
      > You just need to specify your new working directory
  
      ? Enter your rMATS output path (e.g. /User/usr/Tool/SpliceDecoder-main/toy_data)
      > You just need to specify the rMATS output path

      ? Enter your target gene list (e.g. /User/usr/Tool/SpliceDecoder-main/target_genes.tsv)
      > You just need to provide interesting gene list, or enter 'all' if you don’t have one
      > SpliceDeocder will only consider there genes

      ? Enter your GTF file that you used in rMATS with its full path (e.g. /User/usr/Tool/SpliceDecoder-main/toy_data/toy.gtf or /User/usr/Tool/SpliceDecoder-main/toy_data/*.gtf)
      > You just need to specify the full path + GTFfile

      ? Do you want to calculate the effect score? [yes/no]
      > Simply type yes or no. If you type "yes", SpliceDecoder will ask TPM matrix or bamfile path to calculate the effect score

      ? Enter your TPM matrix with full path (e.g. /User/usr/Tool/SpliceDecoder-main/toy_data/tpm.tsv or N)
      > Specify the full path to your TPM matrix, or enter 'N' if you don’t have one
  
      ? Enter your bamlist which should contains bamfile with their full path in each line (e.g. /User/usr/Tool/SpliceDecoder-main/toy_data/bam_list.txt or N)
      > If you don’t have a TPM matrix, specify the full path to your BAM list file, or enter 'N'
  
      ? Enter a species of your data (e.g. human or mouse)
      > You just need to specify the species of your data

      ? Enter a type of GTF (e.g., SR (GENCODE GTF) or LR (Custom GTF) )
      > You just need to specify the type of your GTF

      ? Specify a NMD definition method (e.g., default (55rule) or advanced) )
      > You just need to select one either 'default' or 'advanced'
  
      ? Enter a FDR cut off for your rMATS (float [0-1], default 0.05)
      > Specify rMATS FDR cut off
  
      ? Enter a |dPSI| cut off for your rMATS (float [0-1], default 0.1)
      > Specify rMATS FDR cut off
  
      ? Enter a number of cpu in spliceDecoder job (int [0-?])
      > Specify a number of cpu will be used in your job

* You can reuse a pre-existing config file by copying it:

      cp ${existing_config} project2.config

* Then, update the following fields in the new config: `input`, `Your_GTF`, and `Your_rMATS`

      
  

<br>


## Outputs
```
├── table/
│   ├── *_w_Pfam.txt: Assigned domain information of simulated transcripts (Sim-TXs)
│   └── *_wo_Pfam.txt: Assigned domain information of reference transcripts (Ref-TXs)
├── result/
│   ├── *Main_table.tsv: description
│   ├── *Domain_alt.tsv: description
│   ├── *NMD.tsv: description
│   └── Effect_score.tsv: description
├── figure/
│   ├── mapping_rate.pdf: Mapping rates for each splicing type
│   ├── mat_tx_numbers.pdf: Distribution of Ref-TX for each splicing type
│   ├── splicing_categories_stacked_plot.pdf: description
│   ├── merged_stacked_plot.pdf: description
│   ├── Summary.html: HTML file to make summary pages (pdf_1_page_1.png, pdf_2_page_2.png, pdf_3_page_3.png, and pdf_4_page_4.png)
│   └── consequence: Output directory of visualization script
├── AF2/: Contains AlphaFold2 input (amino acid FASTA)
├── temp/: Contains all intermediate files
├── post_input/: Contains files used in downstream analyiss e.g., visualization and 3D structure generation
├── mapping.stats: Mappeing rates for each splicing type
└── SD.log: The log file
```
* Example of summary HTML
![image](https://github.com/user-attachments/assets/fdfff5a3-b923-45d7-b5aa-4dfcd932c767)

  
<br>

## Details of Outputs
*Example of the `Effect_score.tsv`*
<img width="1040" height="62" alt="image" src="https://github.com/user-attachments/assets/22dcc503-fb7e-4f89-8bfa-0d543e949c63" />

### Key Metrics

- **`LongID`**: DS event ID
- **`Gene symbol`**: Gene symbol
- **`Reference_transcript`**: Matched Transcript (==Ref_TX)
- **`Simulated_event`**: Simulated event (ES = Exon skipping, EI = Exon inclusion, SI = Skipped intron, RI = Retained intron, Can A3/5SS = canonical 3/5' splice site, Alt A3/5SS = alternative 3/5' splice site)
- **`Effect_Score`**: A score to prioritize your DS events [0,2]
- **`Domain_change_rate`**: Average rate of domain changes in Sim-TX compared to Ref-TX [0,1]
- **`Probability_of_NMD`**: NMD **(-1)**, PTC removal **(1)**, No NMD related event **(0)**
- **`Functional_class`**: It contains the following functional classes: GoD (Gain of Domain), LoD (Loss of Domain), NMD, CDS_alts, and UTR_alts
  
### Supplementary Metrics
- `Delta_PSI`: PSI difference (group2 - group1) [-1,1]
- `Transcript_usage`: Proportion of expression of reference transcript for each gene [0,1]
- `ORF`: Used ORF (This file only contains pORF1 which has the highest coding potential)
- `AUG (Ref-Sim)`: Start codon position on the Ref TX and Sim TX (Ref-Sim)
- `Stop_codon (Ref-Sim)`: Stop codon position on the Ref TX and Sim TX (Ref-Sim)
- `Nucleotide_difference`: Coding sequence length difference (Ref TX - Sim TX)
- `5'UTR_difference`: 5' UTR length difference (Ref TX - Sim TX)
- `3'UTR_difference`: 3' UTR length difference (Ref TX - Sim TX)
- `Domain_integrity`: (Sim_domain_length / Ref_domain_length) * 100 [0,inf]
- `Length_of_simulated_tx_domain`: Total domain length of Sim TX
- `Length_of_referece_tx_domain`: Total domain length of Ref TX
- `rMATS_FDR(-log10)`: -Log10 scale FDR, it came from rMATS

---
*Example of the `Domain_alt.tsv`*
![image](https://github.com/user-attachments/assets/266bcdec-d4aa-4e5b-a0bb-3bf0db247a14)

### Key Metrics
- **`DS-TX pair ID`**: It contains, in order Long_ID, Ref-TX ID, and simulated event type
- **`ORF priority`**: A priority of the used reading frame in simulation
- **`Domain information`**: A name of altered domain by the simulated alternative splicing event
- **`Functional_change_ratio (∆L)`**: A difference of functional change ratio for simulated alternative splicing
- **`Change direction`**: It indicates whether the altered domain is a gain (1) or a loss of domain (-1)


---
*Example of the `NMD.tsv`*
<img width="1061" height="177" alt="image" src="https://github.com/user-attachments/assets/c14ca640-4ba6-4feb-a5cd-bef2be03bde2" />

### Key Metrics
- **`LongID`**: Contains, in order Long_ID, Ref-TX ID, and simulated event type
- **`AUG`**: A relative position of AUG on the given transcript (Ref or Sim)
- **`pORF`**: A priority of the used reading frame in simulation
- **`distance(last_exon_junction-stop)`**: Distance between last exon-exon junction and stop codon (calculated by last_exon_junction - stop)
- **`total_domain_length`**: Total domain length of the given transcript (Ref of Sim)
- **`key(Ref/Sim)`**: A type of transcript (Ref of Sim)
- **`NMD_possibility`**: Indicates the possibility of NMD. In default mode, values are HIGH (55nt) or No. In advanced mode, values can be HIGH (55nt), INTERMEDIATE (Long-exon), INTERMEDIATE (Start-proximal), LOW (less 55nt), or No. Only events tagged as HIGH are considered NMD-associated events
- **`contain_PTC`**: Indicates whether the given transcript contains PTC (Y) or not (N)

<br>

## Visualize your alternative splicing simulation
* Based on your Main_table file, you can pcik ceratin DS event to visualize it using this code

      conda activate spliceDecoder
      python code/02-3_v3_Draw_consequence.py \
             --input ${working directory} \
             --splicing_event RI \
             --gene MYLK2 \
             --sim_splicing_event RI \
             --transcript ENSMUST00000028970.7
      python code/02-3_v3_Draw_consequence.py -h  # You can get more details
![RI;ENSMUSG00000027470 9;chr2;+;152919325;152919453;152919454;152920285;152920286;152920438_ENSMUST00000028970 7_splicing_map](https://github.com/user-attachments/assets/9a9e69b1-5ae7-4efb-9228-47c829a0ff40)


* If you want to `remove some information` in figure space, using `ri` option (all categories should be separated by space)

      
      python code/02-3_v3_Draw_consequence.py \
             --input ${working directory} \
             --splicing_event A3SS \
             --gene MYLK2 \
             --sim_splicing_event Can_A3SS \
             --transcript ENSMUST00000195957.4 \
             -ri proteome chain
![image](https://github.com/user-attachments/assets/3306e051-8fa7-47c3-ab25-13c6df061da1)



* All figures will be saved at `${input}/figure/consequence/`

<br>

## Create a 3D Protein structure based on simulated 
* You can use `Make_aa_fa.py` to extract amino acid sequences from your interesting targets
* This function requires the `Effect_score.tsv`, Toy data is not eligible for this function
* You can find the `${input}` and `${Main}` in your `.config` file

      conda activate spliceDecoder
      python code/Make_aa_fa.py \
             -i ${input} \
             -r human \
             -t ENST00000438015.6 \
             -e ES \
             -d ${Main}

* You can copy and paste the amino acid sequences to the Alphafold server (https://alphafoldserver.com) as input

## Please cite this article if you use SpliceDecoder in your research
* https://doi.org/10.1101/2025.10.01.679902
