# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 19:35:47 2020

@author: Daniel Moser
"""

# =============================================================================
# MODULES
# =============================================================================
import os
import re
import math
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO, pairwise2
import matplotlib.pyplot as plt
import logomaker
import Bio.Align.substitution_matrices as matlist
pd.options.mode.chained_assignment = None  # default='warn'

# =============================================================================
# FILES
# =============================================================================
main = r"K:/Sciebo_Koeln/14_Blasts/propeptides_new/"
files = {
    "psiBLAST":
        main+r"PsiBlast_ensemble_allpropeptides_new_E0_5_iter5_max1000000_rep2.csv",
    "Ensemble_Accession_to_Org":
        r"K:/Sciebo_Köln/upload/Ensemble_accession_to_organism.csv",
    "Fasta": r"D:/Bioinformatic/blastdb/ensemble/zipped/ensemble2.fasta",
    "BLASTp": main+r"Blastp_ensemble_allpropeptides_new_E1_max200000.csv",
    "propeptide_to_accession": r"D:/Bioinformatic/sequences/propeptides/new/",
    "summary": main+r"summary/",
    "clustal": main+r"clustal/",
    "clustal2": main+r"clustal_filtered/",
    "clustal3": main+r"clustal_filtered2/",
    "clustal_org": main+r"clustal_filtered_org/"
    }

# =============================================================================
# FUNCTIONS
# =============================================================================

def create_map():
    #matches ensemble accessions to organism
    MAP = pd.read_csv(files["Ensemble_Accession_to_Org"],
                      sep='\t'
                      )
    propeptide_to_accession = {}
    for filename in os.listdir(files["propeptide_to_accession"]):
        if filename.endswith(".fasta"):
            name = filename.replace(".fasta", "")
            name = name.split(" ")[0]
            for record in SeqIO.parse(files["propeptide_to_accession"]+filename,
                                      "fasta"):
                propeptide_to_accession[record.id] = name
    return(MAP, propeptide_to_accession)


def renew_db(BLAST):
    #renews DATABASE from total results
    DATABASE = {}
    PROPEPTIDES = list(BLAST["propeptide"].unique())
    for PROPEPTIDE in PROPEPTIDES:
        TEMP = BLAST[BLAST["propeptide"] == PROPEPTIDE]
        DATABASE[PROPEPTIDE] = TEMP
    return DATABASE


def load_backup(DIR):
    #load backup files
    BLAST = pd.read_excel(DIR, index_col=0)
    DATABASE = renew_db(BLAST)
    return (BLAST, DATABASE)


def process_data(FILE, FILTER):
    #processes BLAST results
    print(FILE)
    BLAST = pd.read_csv(FILE,
                        sep="\t",
                        header=None,
                        comment='#',
                        names=["query id",
                                "subject id",
                                "% identity",
                                "% positives",
                                "% query coverage per subject",
                                "evalue",
                                "bit score",
                                "mismatches",
                                "gaps",
                                "alignment length",
                                "q. start",
                                "q. end",
                                "s. start",
                                "s. end",
                                "subject title"])
    BLAST = BLAST[["query id",
                   "subject id",
                   "% identity",
                   "evalue",
                   "bit score",
                   "% query coverage per subject"
                   ]]
    BLAST.dropna(inplace=True)
    BLAST = BLAST[BLAST["query id"] != "KRH01053"] #subtilase removed from BLAST
    BLAST = BLAST[BLAST['evalue'] <= 0.05] #evalue threshold reduced to 0.05
    BLAST = BLAST.merge(MAP, left_on='subject id', right_on='ACCESSION') #adds organism info to each hit
    del BLAST["ACCESSION"]
    LIST = []
    for ID in BLAST["query id"]:
        LIST.append(propeptide_to_accession[ID]) #adds queried propeptide name to each hit
    BLAST["propeptide"] = LIST    

    print("#hits at start: " + str(len(BLAST)))

    info = {} #collect number of hits at each step
    PROPEPTIDES = sorted(list(BLAST["propeptide"].unique()))  
    for PROPEPTIDE in PROPEPTIDES:
        info[PROPEPTIDE] = PROPEPTIDE + "\t" + str(len(BLAST[BLAST["propeptide"]==PROPEPTIDE]["subject id"]))
   
    #initial quality filter
    BLAST = BLAST[BLAST['% identity'] > 10]
    BLAST = BLAST[BLAST['% query coverage per subject'] > 25]
    print("#hits after quality filter: " + str(len(BLAST)))
    for PROPEPTIDE in PROPEPTIDES:
        info[PROPEPTIDE] += "\t" + str(len(BLAST[BLAST["propeptide"]==PROPEPTIDE]["subject id"]))
    
      
    DATABASE = {}
    for PROPEPTIDE in PROPEPTIDES:
        TEMP = BLAST[BLAST["propeptide"] == PROPEPTIDE]
        TEMP.sort_values(by=['bit score'], inplace=True, ascending=False)
        if FILTER:
            #for psiBLAST: drop reduntant subjects and allow only best hit per organism
            TEMP.drop_duplicates(subset="subject id",
                                 keep='first',
                                 inplace=True)
            TEMP.drop_duplicates(subset="ORG0",
                                 keep='first',
                                 inplace=True,
                                 ignore_index=True)
        IDs = list(TEMP["subject id"])
        new_IDs = []
        for ID in IDs:
            new_IDs.append(ID+"$"+PROPEPTIDE)
        TEMP["unique id"] = new_IDs
        DATABASE[PROPEPTIDE] = TEMP
    LIST = []
    for ID in DATABASE:
        LIST.append(DATABASE[ID])
    BLAST = pd.concat(LIST)
    for PROPEPTIDE in PROPEPTIDES:
        info[PROPEPTIDE] += "\t" + str(len(BLAST[BLAST["propeptide"]==PROPEPTIDE]["subject id"]))
    if FILTER:
        print("#non-redundant hits:" + str(len(BLAST)))
        FASTA = {}
        LEN = {}
        IDS = BLAST["subject id"].to_list()
        for seq_record in SeqIO.parse(files["Fasta"], "fasta"):
            if seq_record.id in IDS:
                FASTA[seq_record.id] = seq_record.seq
                LEN[seq_record.id] = len(FASTA[seq_record.id])
        SEQS = []
        Length = []
        for ID in IDS:
            SEQS.append(str(FASTA[ID]))
            Length.append(LEN[ID])
        BLAST["sequence"] = SEQS
        BLAST["length"] = Length

    DATABASE = renew_db(BLAST) #converts BLAST table to database
    print("\n\n")
    for i in info:
        print(info[i])
    return(BLAST, DATABASE)


def renew_complete(DATABASE):
    # converts DATABASE to total result table
    LIST = []
    for i in DATABASE:
        LIST.append(DATABASE[i])
    BLAST = pd.concat(LIST)
    return BLAST


def safe_data(DATA, NAME):
    #saves FASTA of results
    DATABASE = DATA[1]
    PROPEPTIDES = sorted(DATA[1].keys())
    DIR = files["summary"]+NAME+"/"
    print(DIR)
    for PROPEPTIDE in PROPEPTIDES:
        print(PROPEPTIDE)
        TEMP = DATABASE[PROPEPTIDE]
        texts = ""
        TEMP.set_index("unique id", inplace=True, drop=False)
        for ID in TEMP.index:
            seq = TEMP.loc[ID].at["sequence"]
            text = ">"+str(ID)+"\n" + str(seq)+"\n"
            texts += text
        f = open(DIR+NAME+"_ensemble_propeptides_hit_"+str(PROPEPTIDE)+".fasta", "w")
        f.write(texts)
        f.close()


def read_clustal(DATA, NAME, FILEDIR):
    #reads clustal results
    BLAST = DATA[0]
    clustal = {}
    DIR = FILEDIR+NAME+"/"
    for filename in os.listdir(DIR):
        if filename.endswith(".clustal_num") or filename.endswith(".clustal"):
            print(filename)
            align = AlignIO.read(DIR+filename, "clustal")
            for record in align._records:
                clustal[record.id] = str(record.seq)
    clustal_data = []
    for ID in BLAST["unique id"]:
        if ID in clustal:
            clustal_data.append(clustal[ID])
        else:
            clustal_data.append("")
    BLAST["clustal"] = clustal_data
    DATABASE = renew_db(BLAST)
    return(BLAST, DATABASE, clustal)


def searchable_motif(PROPEPTIDE):
    #creates searchable motif for motif scoring
    motif = peptide_sequences["Peptide"][peptide_sequences["Name"] == PROPEPTIDE].to_string(index=False)
    motif = motif.replace(" ", "")
    motif_s = ""
    for i in motif:
        i = i + "-*"
        motif_s += i
    motif = motif_s[:-2]
    return motif


def clustal_analysis(DATA):
    #analysis of clustal results
    extend = -1
    matrix = matlist.load("BLOSUM62")
    BLAST = DATA[0]
    DATABASE = renew_db(BLAST)
    clustal = DATA[2]
    PROPEPTIDES = list(BLAST["propeptide"].unique())
    HITS = {}
    for PROPEPTIDE in PROPEPTIDES:
        DATABASE[PROPEPTIDE].index = range(0, len(DATABASE[PROPEPTIDE]))
        print(PROPEPTIDE)
        motif = searchable_motif(PROPEPTIDE) #loads motif
        match = re.search(motif, DATABASE[PROPEPTIDE]["clustal"][0]) #searches for motif in clustal
        if match is not None:
            start = int(match.span()[0])
            end = int(match.span()[1])
            matches = {}
            for ID in DATABASE[PROPEPTIDE]["unique id"]:
                matches[ID] = clustal[ID][start:end]
                match0 = match[0].replace("-", "")
                gap = -8  # http://www.genebee.msu.su/blast/blast_help.html
                score0 = pairwise2.align.globalds(match0, match0, matrix, gap, extend)[0].score #alignment and scoring of match
                matchx = matches[ID].replace("-", "")
                if len(matchx) > 0:
                    scorex = pairwise2.align.globalds(match0, matchx, matrix, gap, extend)[0].score
                    score = 100*scorex/score0
                else: 
                    score = 0
                HITS[ID] = [clustal[ID][start:end], matchx, score]
    motifs = []
    motifs_new = []
    motifs_scores = []
    for ID in BLAST["unique id"]:
        if ID in HITS:
            motifs.append(HITS[ID][0])
            motifs_new.append(HITS[ID][1])
            motifs_scores.append(HITS[ID][2])
        else:
            motifs.append("")
            motifs_new.append("")
            motifs_scores.append(0)
    BLAST["motif"] = motifs
    BLAST["motif_new"] = motifs_new
    BLAST["motif_score"] = motifs_scores
    DATABASE = renew_db(BLAST)
    return(BLAST, DATABASE, clustal)


def compare(BLASTp, psiBLAST):
    #compares results of BLASTp and psiBLAST
    psiBLAST = psiBLAST[0]
    BLASTp = BLASTp[0]
    IDs_BLASTp = list(BLASTp["unique id"])
    IDs_psiBLAST = list(psiBLAST["unique id"])
    result_psiBLAST = []
    # check if BLASTp in psiBLAST
    for ID in IDs_psiBLAST:
        if ID in IDs_BLASTp:
            result_psiBLAST.append(True)
        else:
            result_psiBLAST.append(False)
    psiBLAST["BLASTp"] = result_psiBLAST
    psiBLAST_db = renew_db(psiBLAST)
    BLASTp_db = renew_db(BLASTp)
    return([psiBLAST, psiBLAST_db], [BLASTp, BLASTp_db])


def filter_data(DATA, threshold):
    #filters DATA by motif score
    BLAST = DATA[0][DATA[0]["motif_score"] > threshold]
    DATABASE = renew_db(BLAST)
    return (BLAST, DATABASE)


def overview_id_motif(DATA, selection, folder):
    #creates overview plots of hits for each propeptide
    PROPEPTIDES = selection
    print(PROPEPTIDES)
    x = int(np.ceil(len(PROPEPTIDES)/2))
    fig = plt.figure(figsize=(9, 9))
    gs = fig.add_gridspec(x, 2, hspace=0.1, wspace=0.04)
    axs = gs.subplots(sharex=True, sharey=True)
    print(axs)

    for PROPEPTIDE in PROPEPTIDES:
        print(PROPEPTIDE)
        i = PROPEPTIDES.index(PROPEPTIDE)
        print(i)
        if i % 2 == 0:
            c = 0
        else:
            c = 1
        r = int(np.floor(i/2))
        TEMP = DATA[0][DATA[0]["propeptide"] == PROPEPTIDE]
        a = list(TEMP["% identity"])
        b = list(TEMP["motif_score"])
        axs[r,c].plot(a, b, ".", color="black", markersize=8)
        

        
        TEMPx = TEMP[TEMP["BLASTp"] == True]
        a = TEMPx["% identity"]
        b = TEMPx["motif_score"]
        axs[r,c].plot(a, b, 'x', color="red", alpha=0.5, markersize=8)
        axs[r,c].set_title(PROPEPTIDE, y = 0.85, size=14, fontweight="bold")
        
        TEMPx = TEMP[TEMP["ORG1"] == "Zea"]
        a = TEMPx["% identity"]
        b = TEMPx["motif_score"]
        axs[r,c].plot(a, b, '.', color="lawngreen", markersize=8)
    
    for ax in axs.flat:

        ax.set(yticks=np.arange(0,121,20))
        ax.set_xlabel("% Protein identity", fontsize=12, fontweight="bold")
        ax.set_ylabel("% Motif score", fontsize=12, fontweight="bold")
        ax.label_outer()

def plot_overview(DATA, x, limit, xlimit, folder):
    #creates overview plots of hits for each propeptide, selected column vs. motif score
    PROPEPTIDES = list(DATA[0]["propeptide"].unique())
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 14,
            }
    for PROPEPTIDE in PROPEPTIDES:
        TEMP = DATA[0][DATA[0]["propeptide"] == PROPEPTIDE]
        a = list(TEMP[x])
        b = list(TEMP["motif_score"])
        print(PROPEPTIDE)
        plt.figure(figsize=(3, 3))
        plt.plot(a, b, '.', color="black")

        TEMPx = TEMP[TEMP["ORG1"] == "Zea"]
        a = TEMPx[x]
        b = TEMPx["motif_score"]
        plt.plot(a, b, '.', color="blue")
        
        TEMPx = TEMP[TEMP["BLASTp"] == True]
        a = TEMPx[x]
        b = TEMPx["motif_score"]
        plt.plot(a, b, 'x', color="red", alpha=0.5)
        plt.title(PROPEPTIDE, fontdict=font, y=0.05, x=0.5)
        plt.ylim(0, 110)
        if xlimit is True:
            plt.xlim(0, 110)
        plt.xticks([])
        plt.ylabel("Motif score [%]", fontdict=font)
        plt.tight_layout()
        plt.savefig(main+
                    folder+'/'+PROPEPTIDE+"motif-score_vs_"+x+".png", dpi=600)
        plt.show()


def create_fasta(DATA):
    PROPEPTIDES = list(DATA[0]["propeptide"].unique())
    DB = DATA[1]
    FASTAS = {}
    for PROPEPTIDE in PROPEPTIDES:
        TEMP = DB[PROPEPTIDE]
        IDS = list(TEMP["unique id"])
        TEMP = TEMP.set_index("unique id")
        FASTA = []
        for ID in IDS:
            motif = TEMP.loc[ID].at["motif"]
            FASTA.append(motif)
        FASTAS[PROPEPTIDE] = FASTA
    return FASTAS


def create_weblogo(FASTAS, folder):
    #creates motif logos
    PROPEPTIDES = list(FASTAS.keys())
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 16,
            }
    for PROPEPTIDE in PROPEPTIDES:
        print(PROPEPTIDE)
        if len(FASTAS[PROPEPTIDE]) <= 1:
            print(PROPEPTIDE+" has only 1 or less sequences")
        else:
            df = {}
            df = logomaker.alignment_to_matrix(FASTAS[PROPEPTIDE], to_type="counts")
            LOGO = logomaker.Logo(df,
                                  color_scheme="dmslogo_funcgroup", #or chemistry
                                  font_name="Times New Roman",
                                  stack_order="big_on_top",
                                  show_spines=True,
                                  vpad=0.05)
            LOGO.style_xticks(rotation=0, fmt='%d', anchor=0, spacing=5)
            LOGO.draw(clear=False)
            plt.title(PROPEPTIDE, fontdict=font)
            plt.xticks(size=14)
            plt.yticks(size=14)
            plt.xlabel("Position", fontdict=font)
            plt.ylabel("Count", fontdict=font)
            plt.tight_layout()
            plt.savefig(main+
                        folder+'/'+PROPEPTIDE+".png", dpi=600)
            plt.show()


def prepare_clustal(DATA, folder):
    #prepares FASTA for clustal
    DATABASE = DATA[1]
    DIR = main+folder+"/"
    for PROPEPTIDE in DATA[1]:
        texts = ""
        for i in DATABASE[PROPEPTIDE]["unique id"]:
            seq = DATABASE[PROPEPTIDE][DATABASE[PROPEPTIDE]["unique id"] == i]["sequence"].values[0]
            text = ">" + str(i) + "\n" + str(seq) + "\n"
            texts += text
        f = open(DIR+"filtered_"+str(PROPEPTIDE)+".fasta", "w")
        f.write(texts)
        f.close()


def plot_overview2(DATA):
    #other type of overview plot, length vs identity
    DATABASE = DATA[1]
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 12,
            }
    x = "length"
    y = "% identity"
    
    for PROPEPTIDE in DATABASE:
        plt.figure(figsize=(2,2))
        TEMP = DATABASE[PROPEPTIDE]
        b = TEMP[y]
        a = TEMP[x]
        a_max = math.ceil(max(a)*1.1)
        plt.xticks(ticks=np.arange(0, a_max, a_max/5), rotation=90)
        plt.plot(a, b, '.', color="black")
        TEMPx = TEMP[TEMP["length_criteria"] == False]
        a = TEMPx[x]
        b = TEMPx[y]
        plt.plot(a, b, '.', color="red")
        plt.title(PROPEPTIDE, fontdict=font)
        plt.ylim(0, 110)
        
        #plt.xlim(0)
        plt.xlabel("Length", fontdict=font)
        plt.ylabel("% Identity", fontdict=font)
        plt.tight_layout()
        
        plt.savefig(main+'plots/length_selection/length_selection_' +
                        PROPEPTIDE+".png", dpi=600)
        plt.show()


def rel_length(DATA):
    #calculates relative length of each hit
    #compared to query
    # and filters
    DATABASE = DATA[1]
    length_criteria = {}
    x = []
    for PROPEPTIDE in DATABASE:
        TEMP = DATABASE[PROPEPTIDE]
        TEMP.set_index("unique id", inplace=True, drop=False)
        LENGTH75 = TEMP["length"][TEMP["% identity"] > 75]
        minima = min(LENGTH75)*0.5
        maxima = max(LENGTH75)*1.5
        for ID in TEMP.index:
            LENGTH_i = TEMP.loc[ID].at["length"]
            if LENGTH_i < maxima and LENGTH_i > minima:
                length_criteria[ID] = True
            else:
                length_criteria[ID] = False
    for ID in DATA[0]["unique id"]:
        x.append(length_criteria[ID])
    DATA[0]["length_criteria"] = x
    BLAST = DATA[0][DATA[0]["length_criteria"] == True]
    DATABASE = renew_db(BLAST)
    for PROPEPTIDE in sorted(list(DATABASE.keys())):
        print(PROPEPTIDE + "\t" + str(len(DATABASE[PROPEPTIDE]["subject id"])))
    return(BLAST, DATABASE)


def create_clustal_cmd(DATA, NAME):
    #performes clustal
    from Bio.Align.Applications import ClustalOmegaCommandline
    #if NAME  == "clustal":
    DIR_in = main+"summary/psiBLAST/"
    #else:
    #    DIR_in = main+"for_clustal_filtered/"
    DIR_out = main+NAME+"/psiBLAST/"
    command = "echo START of ClustalO"
    PROPEPTIDES = DATA[1].keys()
    for PROPEPTIDE in PROPEPTIDES:
        print(PROPEPTIDE)
        #if NAME == "clustal":
        in_file = DIR_in+"psiBLAST_ensemble_propeptides_hit_" + PROPEPTIDE + ".fasta"
        #else:
            #in_file = DIR_in + "filtered_"+PROPEPTIDE+".fasta"
        out_file = DIR_out + PROPEPTIDE + ".clustal"
        clustalomega_cline = ClustalOmegaCommandline(
            cmd=r"K:\clustal-omega\clustal-omega-1.2.2-win64\clustalo",
            infile=in_file,
            outfile=out_file,
            outfmt="clustal",
            verbose=True,
            outputorder="input-order",
            auto=True)
        command += " & echo " + PROPEPTIDE
        command += " & " + str(clustalomega_cline) + " --force & echo ."
    command += " & echo END of ClustalO"
    print(command)
    os.system('cmd /k "'+command+'"')


def histogram(DATA, PROPEPTIDE):
    #creates histogram of hits after each step
    if PROPEPTIDE is False:
        if type(DATA[0]) == tuple:
            for i in DATA:
                x = i[0]["% identity"]
                plt.hist(x, bins=100)
        else:
            x = DATA[0]["% identity"]
            plt.hist(x, bins=100)
        plt.ylabel("Number of hits", fontweight='bold')
        plt.xlabel("% Identity", fontweight='bold')
        plt.xticks(np.arange(0, 101, 10))
        plt.legend(["Pre-filtered hits", "After filtering by length", "After iterative motif filtering"])
        plt.tight_layout()
        plt.savefig(main+'plots/histogram.png', dpi=600)
        plt.show()
    else:
        for PROPEPTIDE in DATA[1].keys():
            x = DATA[1][PROPEPTIDE]["% identity"]
            plt.hist(x, bins=100)
            plt.ylabel("Number of hits", fontweight='bold')
            plt.xlabel("% Identity", fontweight='bold')
            plt.xticks(np.arange(0, 101, 10))
            plt.title(PROPEPTIDE)
            plt.tight_layout()
            plt.savefig(main+'plots/histogram_' +
                        PROPEPTIDE+".png", dpi=600)
            plt.show()


def violinplot(DATA, key, zoom, name):
    #violin plot of different factors
    TEMP = DATA[1]
    x = []
    label_propeptides = []
    PROPEPTIDES = sorted(TEMP.keys())
    for PROPEPTIDE in PROPEPTIDES:
        x.append(TEMP[PROPEPTIDE][key])
        label_propeptides.append(PROPEPTIDE + " ("+str(len(TEMP[PROPEPTIDE][key]))+")")
    fig, ax = plt.subplots()
    plot = ax.violinplot(x,
                         showmedians=True,
                         showextrema=False,
                         vert=False,
                         widths=0.7,
                         points=500)
    
    for pc in plot['bodies']:
        pc.set_facecolor('grey')
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
    plot['cmedians'].set_edgecolor("red")
    
    plt.title("")
    ax.set_title('')
    if key == "motif_score":
        ax.set_xticks(range(0, 101, 10))
        plt.xticks(rotation='vertical')
        plt.xlabel("Motif score [%]", fontweight='bold')
    if key == "length":
        plt.xticks(rotation='vertical')
        if zoom == False:
            plt.xticks(range(0, 6001, 500))        
        plt.xlabel("Protein length", fontweight='bold')
    if key == "% identity":
        plt.xticks(range(0, 101, 10), rotation="vertical")
        plt.xlabel("Protein identity [%]", fontweight='bold')
    ax.set_yticks(range(1, len(PROPEPTIDES)+1))
    ax.set_yticklabels(label_propeptides)
    plt.tight_layout()    
    plt.savefig(main+'plots/violinplot_' +
                str(name)+"_"+str(key)+"_zoom_"+str(zoom)+".png", dpi=600)
    plt.show()
    
def violinplot2(DATA):
    #violin plot of motif score before/after length filtering
    PROPEPTIDES = list(DATA[0][1].keys())
    PROPEPTIDES.remove("PreproHypSys")
    plt.figure(figsize=[10,10], dpi=600)
    fig, axs = plt.subplots(len(PROPEPTIDES), 1, sharex=True)
    
    for PROPEPTIDE in PROPEPTIDES:
        x = []
        for DATAi in DATA:
            x.append(list(DATAi[1][PROPEPTIDE]["motif_score"]))
        i_p = PROPEPTIDES.index(PROPEPTIDE)
        i = 0
        axs[i_p].violinplot(x[i],
                          showmedians=True,
                          showextrema=False,
                          vert=False,
                          widths=0.5,
                          points=100)
        axs[i_p].violinplot(x[i+1],
                          showmedians=True,
                          showextrema=False,
                          vert=False,
                          widths=0.5,
                          points=100)
        axs[i_p].set_yticks([])
        x0 = str(len(x[i]))
        x1 = str(len(x[i+1]))
        title = PROPEPTIDE + " (" + x0 + "|" + x1 + ")"
        axs[i_p].set_title(title, x= -0.12, y=0.75, pad=-14, fontweight="bold")
    
    fig.set_size_inches(10, 10)
    fig.legend(axs,     
            labels=["Distribution\nbefore filtering", 
                    "Median\nbefore filtering",
                    "Distribution\nafter filtering",
                    "Median\nafter filtering"],
            bbox_to_anchor=(0., 0.94, 1., .102),
            #loc="lower right", 
            #borderaxespad=0.1,
            ncol=4, mode="expand"
            )
    plt.xlabel("Motif score [%]", fontweight="bold", fontsize="medium")
    
    plt.tight_layout()
    plt.savefig(main+'plots/violinplot_motifscoreBeforeAfterMotifScoring.png',
                dpi=600,
                format="png",
                bbox_inches='tight')
    


def motif_hitlist(DATA):
    #creates a list of hits per propeptide
    DATABASE = DATA[1]
    for PROPEPTIDE in DATABASE.keys():
        print(PROPEPTIDE)
        print("Alignment\t% motif score\t%identity\t organism")
        for i in list(DATABASE[PROPEPTIDE].index):
            print(str(DATABASE[PROPEPTIDE].loc[i].at["motif"])+"\t"
                  + str(DATABASE[PROPEPTIDE].loc[i].at["motif_score"].round(2))+"\t"
                  + str(DATABASE[PROPEPTIDE].loc[i].at["% identity"].round(2))+"\t"
                  + str(DATABASE[PROPEPTIDE].loc[i].at["ORG0"])+"\t")
        print("\n\n")


def Zea_hits(DATA):
    #creates a list of all hits found in Z. mays
    DATABASE = DATA[1]
    print("Propeptide\tOrganism\tID\tPeptide hormone motif\tIdentity [%]\tMotif score [%]\tLength\tHit in BLASTp?")
    for PROPEPTIDE in DATABASE.keys():
        if list(DATABASE[PROPEPTIDE]["ORG0"]).count("Zea_mays") >= 1 and PROPEPTIDE != "ZmProPEP1":           
            list0 = list(DATABASE[PROPEPTIDE].index)
            i0 = list0[0]
            for i in list(DATABASE[PROPEPTIDE].index):
                ORG = DATABASE[PROPEPTIDE].loc[i].at["ORG0"]
                if ORG == "Zea_mays" or ORG == "Arabidopsis_thaliana" or i == i0:
                    print(str(DATABASE[PROPEPTIDE].loc[i].at["propeptide"]) + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["ORG0"]).replace("_"," ") + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["subject id"]) + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["motif"]) + "\t" #"-" need to me removed by hand
                          + str(DATABASE[PROPEPTIDE].loc[i].at["% identity"].round(2)) + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["motif_score"].round(2)) + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["length"]) + "\t"
                          + str(DATABASE[PROPEPTIDE].loc[i].at["BLASTp"]) + "\t"
                          )


def info(DATA):
    print("total data\t" + str(len(DATA[0])))
    for PROPEPTIDE in DATA[1]:
        print(PROPEPTIDE+"\t"+str(len(DATA[1][PROPEPTIDE])))


def info_orginaldata(FILE, FILTER):
    #not used anymore
    print(FILE+"\n")
    BLAST = pd.read_csv(FILE,
                        sep="\t",
                        header=None,
                        comment='#',
                        names=["query id",
                                "subject id",
                                "% identity",
                                "% positives",
                                "% query coverage per subject",
                                "evalue",
                                "bit score",
                                "mismatches",
                                "gaps",
                                "alignment length",
                                "q. start",
                                "q. end",
                                "s. start",
                                "s. end",
                                "subject title"])
    BLAST = BLAST[["query id",
                   "subject id",
                   "% identity",
                   "evalue",
                   "bit score",
                   "% query coverage per subject"
                   ]]
    BLAST.dropna(inplace=True)
    BLAST = BLAST.merge(MAP, left_on='subject id', right_on='ACCESSION')
    del BLAST["ACCESSION"]
    LIST = []
    for i in BLAST["query id"]:
        LIST.append(propeptide_to_accession[i])
    BLAST["propeptide"] = LIST
    PROPEPTIDES = list(BLAST["propeptide"].unique())
    DATABASE = {}
    print("total hits\t" + str(len(BLAST)) + "\n")
    if FILTER:
        text = "Propeptide\thits\tw/o redundant\tbest per organism\tfiltered\n"
    else:
        text = "Propeptide\thits\n"
    for PROPEPTIDE in PROPEPTIDES:
        TEMP = BLAST[BLAST["propeptide"] == PROPEPTIDE]
        TEMP.sort_values(by=['bit score'], inplace=True, ascending=False)

        if FILTER:
            text += PROPEPTIDE + "\t" + str(len(TEMP))
            TEMP.drop_duplicates(subset="subject id",
                                 keep='first',
                                 inplace=True)
            text += "\t" + str(len(TEMP))
            TEMP.drop_duplicates(subset="ORG0",
                                 keep='first',
                                 inplace=True,
                                 ignore_index=True)
            text += "\t" + str(len(TEMP))
            TEMP = TEMP[TEMP['% identity'] > 10]
            TEMP = TEMP[TEMP['% query coverage per subject'] > 25]
            text += "\t" + str(len(TEMP)) + "\n"
        else:
            text += PROPEPTIDE + "\t" + str(len(TEMP)) + "\n"
        IDs = list(TEMP["subject id"])
        new_IDs = []
        for ID in IDs:
            new_IDs.append(ID+"$"+PROPEPTIDE)
        TEMP["unique id"] = new_IDs
        DATABASE[PROPEPTIDE] = TEMP
    LIST = []
    for i in DATABASE:
        LIST.append(DATABASE[i])
    BLAST = pd.concat(LIST)
    text += "\nfiltered hits\t" + str(len(BLAST)) + "\n"
    print(text)


def collect_infos():
    ORGANISMS = sorted(list(psiBLAST_filtered3[0]["ORG0"].unique()))
    print("#proteins in BLAST database:\t" + str(len(MAP)))
    print("Organisms in BLAST database:\t" + str(len(ORGANISMS)))
    for ORGANISM in ORGANISMS:
        print(ORGANISM.replace("_", " "))
    return(ORGANISMS)

def filter_organisms(DATA):
    #filter out genera duplicates
    DATABASE = DATA[1]
    for PROPEPTIDE in DATABASE:
        TEMP = DATABASE[PROPEPTIDE]
        TEMP.sort_values(by=['motif_score'], inplace=True, ascending=False)
        TEMP.drop_duplicates(subset="ORG1",
                            keep='first',
                            inplace=True,
                            ignore_index=True)
        DATABASE[PROPEPTIDE] = TEMP
    BLAST = renew_complete(DATABASE)
    return(BLAST, DATABASE)

def clustal_evaluation(data_list):
    #evaluation plot of clustal steps
    PROPEPTIDES = sorted(data_list[0][1].keys())
    evaluation_all = []
    for data in data_list:
        evaluation = []
        data = data[1]
        #PROPEPTIDES = sorted(data.keys())
        for PROPEPTIDE in PROPEPTIDES:
            if PROPEPTIDE in data.keys():
                if len(data[PROPEPTIDE]["clustal"])==1:
                    evaluation.append(0)
                else:
                    rel = len(data[PROPEPTIDE]["clustal"][0])/data[PROPEPTIDE]["length"][0]
                    evaluation.append(rel)
            else: 
                evaluation.append(0)
        evaluation_all.append(evaluation) 
    
    labels = PROPEPTIDES
    al1 = evaluation_all[0]
    al2 = evaluation_all[1]
    al3 = evaluation_all[2]
    x = np.arange(len(labels))
    width = 0.2
    fig, ax = plt.subplots()
    ax.bar(x - width, al1, width, label='1st Aligment')
    ax.bar(x , al2, width, label='2nd Aligment')
    ax.bar(x + width, al3, width, label='3rd Aligment')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation="vertical")  
    plt.xlabel("Propeptides", fontsize="large", fontweight='bold')
    plt.ylabel("Relative alignment length", fontsize="large", fontweight='bold')
    plt.legend(labels=["1st Alignment", "2nd Alignment", "3rd Alignment"],
               #bbox_to_anchor=(1.00, 0.75),
               loc='upper right',
               fontsize = "small")
    fig.tight_layout()
    plt.savefig(main+'plots/barplot_alignmentlength.png', dpi=600)
    plt.show()

MAP, propeptide_to_accession = create_map() #loading propeptide to accession table
peptide_sequences = pd.read_excel(main+'propeptide peptide sequences.xlsx')
#processing psiBLAST data + intitial filtering
psiBLAST = process_data(files["psiBLAST"], True)
psiBLAST[0].to_excel(main+"summary/psiBLAST_ensemble_propeptides_hits.xlsx")
psiBLAST = load_backup(main+"summary/psiBLAST_ensemble_propeptides_hits.xlsx")
#processing BLASTp data + intitial filtering
BLASTp = process_data(files["BLASTp"], False)
BLASTp[0].to_excel(main+"summary/BLASTp_ensemble_propeptides_hits.xlsx")
BLASTp = load_backup(main+"summary/BLASTp_ensemble_propeptides_hits.xlsx")

#Length Filtering
psiBLAST1 = rel_length(psiBLAST)
psiBLAST1[0].to_excel(main+"summary/psiBLAST1_ensemble_propeptides_hits.xlsx")
psiBLAST1 = load_backup(main+"summary/psiBLAST1_ensemble_propeptides_hits.xlsx")
safe_data(psiBLAST1, "psiBLAST1") #saving fastas files

#first iterative motif scoring
create_clustal_cmd(psiBLAST1, "clustal")
psiBLAST2 = read_clustal(psiBLAST1, "psiBLAST", files["clustal"])
DIR = main+"summary/psiBLAST2_ensemble_propeptides_hits.xlsx"
psiBLAST2[0].to_excel(DIR)
psiBLAST2 = load_backup(main+"summary/psiBLAST2_ensemble_propeptides_hits.xlsx")
psiBLAST3 = clustal_analysis(psiBLAST2)
DIR = main+"summary/psiBLAST3_ensemble_propeptides_hits.xlsx"
psiBLAST3[0].to_excel(DIR)
psiBLAST3 = load_backup(main+"summary/psiBLAST3_ensemble_propeptides_hits.xlsx")

psiBLAST4, BLASTp1 = compare(BLASTp, psiBLAST3) #compares psiBLAST and BLASTp

DIR = main+"summary/psiBLAST4_ensemble_propeptides_hits.xlsx"
psiBLAST4[0].to_excel(DIR)

psiBLAST4 = load_backup(main+"summary/psiBLAST4_ensemble_propeptides_hits.xlsx")
psiBLAST_filtered = filter_data(psiBLAST4, 10) # filter by motif score

#create logo of motifs
psiBLAST_motif_FASTA = create_fasta(psiBLAST_filtered)
create_weblogo(psiBLAST_motif_FASTA, "logos")

#second iterative motif scoring
prepare_clustal(psiBLAST_filtered, "for_clustal_filtered")
create_clustal_cmd(psiBLAST_filtered, "clustal_filtered")
psiBLAST5 = read_clustal(psiBLAST_filtered, "psiBLAST", files["clustal2"])
DIR = main+"summary/psiBLAST5_ensemble_propeptides_hits.xlsx"
psiBLAST5[0].to_excel(DIR)
psiBLAST6 = clustal_analysis(psiBLAST5)
DIR = main+"summary/psiBLAST6_ensemble_propeptides_hits.xlsx"
psiBLAST6[0].to_excel(DIR)
psiBLAST6 = load_backup(main+"summary/psiBLAST6_ensemble_propeptides_hits.xlsx")
psiBLAST_filtered2 = filter_data(psiBLAST6, 10) #filtering motif score

#third iterative motif scoring
psiBLAST_2_motif_FASTA = create_fasta(psiBLAST_filtered2)
prepare_clustal(psiBLAST_filtered2, "for_clustal_filtered2")
create_clustal_cmd(psiBLAST_filtered2, "clustal_filtered2")
psiBLAST7 = read_clustal(psiBLAST_filtered2, "psiBLAST", files["clustal3"])
psiBLAST8 = clustal_analysis(psiBLAST7)
DIR = main+"summary/psiBLAST8_ensemble_propeptides_hits.xlsx"
psiBLAST8[0].to_excel(DIR)
psiBLAST8 = load_backup(main+"summary/psiBLAST8_ensemble_propeptides_hits.xlsx")
psiBLAST_filtered3 = filter_data(psiBLAST8, 0.10) #filter motif score

DIR = main+"summary/psiBLAST_filtered3_ensemble_propeptides_hits.xlsx"
psiBLAST_filtered3[0].to_excel(DIR)
psiBLAST_filtered3 = load_backup(main+"summary/psiBLAST_filtered3_ensemble_propeptides_hits.xlsx")

psiBLAST_3_motif_FASTA = create_fasta(psiBLAST_filtered3)

create_weblogo(psiBLAST_3_motif_FASTA, "logos_filtered2") #final logo

#gathering information from results
info(psiBLAST_filtered3)
organisms = collect_infos()
info_orginaldata(files["psiBLAST"], True)
info_orginaldata(files["BLASTp"], False)

#creation of plots
#histogram steps
histogram([psiBLAST, psiBLAST1, psiBLAST_filtered3], False)

#Violin plots final
violinplot(psiBLAST_filtered3, "motif_score", False, "motif_score3")
violinplot(psiBLAST_filtered3, "% identity", False, "idenity3")
violinplot(psiBLAST_filtered3, "length", True, "length3")

#Comparison Length filtering
violinplot(psiBLAST, "length", False, "initial")
violinplot(psiBLAST1, "length", False, "length_filtered")

#Identity vs. motif score
plot_overview(psiBLAST_filtered3, "% identity", 10, True, "plots")

#Relative alignment length
clustal_evaluation([psiBLAST_filtered, psiBLAST_filtered2, psiBLAST_filtered3])

#Before / after motif filtering
violinplot2([psiBLAST4, psiBLAST_filtered3])

#list of  Zea mays hits
Zea_hits(psiBLAST_filtered3)

#identity vs motif score of selection combined in one plot
overview_id_motif(psiBLAST_filtered3,
                  ["PSK1",
                   "PrePIP1",
                   "ProRALF23",
                   "PSY1",
                   "IRP",
                   "ProPEP1"],
                  "")



#web logo without redundant organisms in genera
psiBLAST9 = filter_organisms(psiBLAST_filtered3)
prepare_clustal(psiBLAST9, "for_clustal_filtered_org")
create_clustal_cmd(psiBLAST9, "clustal_filtered_org")
psiBLAST10 = read_clustal(psiBLAST9, "psiBLAST", files["clustal_org"])
psiBLAST11 = clustal_analysis(psiBLAST10)
DIR = main+"summary/psiBLAST11_ensemble_propeptides_hits.xlsx"
psiBLAST11[0].to_excel(DIR)
psiBLAST_4_motif_FASTA = create_fasta(psiBLAST11)   
create_weblogo(psiBLAST_4_motif_FASTA, "logos_organisms_filtered")

#final motif hit list
motif_hitlist(psiBLAST_filtered3) 
print("Propeptide\tID\tOrganism\tLength\tPetide hormone")
for PROPEPTIDE in sorted(psiBLAST_filtered3[1].keys()):
    ID = str(psiBLAST_filtered3[1][PROPEPTIDE]["query id"][0])
    ORG0 = str(psiBLAST_filtered3[1][PROPEPTIDE]["ORG0"][0])
    Length = str(psiBLAST_filtered3[1][PROPEPTIDE]["length"][0])
    Peptide = str(psiBLAST_filtered3[1][PROPEPTIDE]["motif_new"][0])
    print(str(PROPEPTIDE) +
          "\t" + ID +
          "\t" + ORG0 +
          "\t" + Length +
          "\t" + Peptide
          )

# plot_overview(psiBLAST_filtered3, "length", 5, False, "plots")

# plot_overview(psiBLAST6, "% identity", 5, True, "plots_filtered")
# plot_overview2(psiBLAST1)
# create_weblogo(psiBLAST_2_motif_FASTA, "logos_filtered")
# histogram(psiBLAST, False)
# plot_overview2(psiBLAST6)
# histogram(psiBLAST8, True)

#another info table
def info2(PROPEPTIDE):
    print("Alignment\t% motif score\t%identity\t organism")
    
    for i in list(psiBLAST4[1][PROPEPTIDE].index):
        print(psiBLAST4[1][PROPEPTIDE].loc[i].at["motif"]+"\t"
              +str(psiBLAST4[1][PROPEPTIDE].loc[i].at["motif_score"].round(2))+"\t"
              +str(psiBLAST4[1][PROPEPTIDE].loc[i].at["% identity"].round(2))+"\t"
              +psiBLAST4[1][PROPEPTIDE].loc[i].at["ORG0"]+"\t")

info2("Subtilase")
info(psiBLAST9)

#longer RALF23 motif test
def motif_longer(DATA, PROPEPTIDE):
    TEMP = DATA[1][PROPEPTIDE]
    motif0 = TEMP["motif"][0]
    clustal0 = TEMP["clustal"][0]
    print(motif0)
    hit0 = re.search(motif0, clustal0)
    start = int(hit0.span()[0])
    end = int(hit0.span()[1])
    new = []
    for i in list(TEMP.index):
        new.append(TEMP.loc[i].at["clustal"][start:end])
    TEMP["new"] = new
    return(TEMP)
Test = filter_organisms(psiBLAST_filtered3)
Test = motif_longer(Test, "ProRALF23")


#BUSCO analysis
#/media/alle/Bioinformatic/Bioinformatic/busco
DIRi = "D:/Bioinformatic/" #/media/alle/Bioinformatic/Bioinformatic/
DIR = DIRi+"busco"
#"D:/Bioinformatic/busco"
folders = os.listdir(DIR)
folders.remove("busco_downloads")
org = [s.split(".")[0] for s in folders]

#organisms selected for phylogenetic tree
organisms = [
     'Arabidopsis_lyrata',
     'Arabidopsis_thaliana',
     'Arabis_alpina',
     'Beta_vulgaris',
     'Brachypodium_distachyon',
     'Brassica_napus',
     'Brassica_rapa',
     'Glycine_max',
     'Gossypium_raimondii',
     'Hordeum_vulgare',
     'Hordeum_vulgare_goldenpromise',
     'Malus_domestica_golden',
     'Nicotiana_attenuata',
     'Oryza_indica',
     'Oryza_sativa',
     'Physcomitrella_patens',
     'Oryza_nivara',
     'Solanum_lycopersicum',
     'Solanum_tuberosum',
     'Sorghum_bicolor',
     'Trifolium_pratense',
     'Triticum_aestivum',
     'Triticum_dicoccoides',
     'Zea_mays'
     ]

results = {}
for folder in folders:
    if folder.split(".")[0] in organisms:
        p = DIR+"/"+folder+"/run_viridiplantae_odb10/full_table.tsv"
        print(folder.split(".")[0])
        table = pd.read_csv(p, sep="\t", skiprows=2)
        results[folder.split(".")[0]] = table
        Stati = list(table["Status"])
        S = Stati.count("Complete")
        D = Stati.count("Duplicated")
        F = Stati.count("Fragmented")
        M = Stati.count("Missing")
        print(S, D, F, M)

#collects all single and multi copy BUSCO gene ids
files = []
for folder in folders:
    if folder.split(".")[0] in organisms:
        p = DIR+"/"+folder+"/run_viridiplantae_odb10/busco_sequences/single_copy_busco_sequences"
        single = os.listdir(p)
        p = DIR+"/"+folder+"/run_viridiplantae_odb10/busco_sequences/multi_copy_busco_sequences"
        duplicates = os.listdir(p)
        genes = single + duplicates 
        fil = {names for names in genes if names.endswith(".faa")}
        files.append(fil)
  
tree_genes = set.intersection(*files) #creates overlap between all organism's BUSCO genes
print(len(tree_genes))

#performs clustal alignment of each BUSCO gene
hits = {}
from Bio.Align.Applications import ClustalOmegaCommandline
for gene in tree_genes:
    #print(gene)
    FASTA = ""
    for folder in folders:
        organism = folder.split(".")[0]
        if organism in organisms:
            pS = DIR+"/"+folder+"/run_viridiplantae_odb10/busco_sequences/single_copy_busco_sequences/"
            pM = DIR+"/"+folder+"/run_viridiplantae_odb10/busco_sequences/multi_copy_busco_sequences/"
            if gene in os.listdir(pS):
                p = pS
            else: p = pM
            p = p + gene
            records = list(SeqIO.parse(p, "fasta"))
            record = records[0]
            FASTA += ">" + organism + "\n"
            FASTA += str(record.seq) + "\n"
            hits[gene] = FASTA
    fasta_dir = "/media/alle/Bioinformatic/Bioinformatic/busco_clustalo/"
    f = open(fasta_dir + gene.replace(".faa", "") + ".fasta", "w")
    f.write(FASTA)
    f.close()
    
command = "echo START of ClustalO"
for gene in tree_genes:
    clustal_DIR = "/media/alle/Bioinformatic/Bioinformatic/busco_clustalo/"
    in_file = clustal_DIR + gene.replace(".faa", "") + ".fasta"
    out_file = clustal_DIR + gene.replace(".faa", "") + ".clustal"
    clustalomega_cline = ClustalOmegaCommandline(
        #cmd=r"K:\clustal-omega\clustal-omega-1.2.2-win64\clustalo",
        infile=in_file,
        outfile=out_file,
        outfmt="clustal",
        verbose=True,
        outputorder="input-order",
        auto=True)
    command += " && " + str(clustalomega_cline) + " --force && echo ."
command += " && echo END of ClustalO"
print(command)

# concatenation of each alignment for phylogenetic tree creation
align_org = {}
DIR = "/media/alle/Bioinformatic/Bioinformatic/busco_clustalo/"
for gene in tree_genes:
    gene = gene.replace(".faa", "")
    for alignment in AlignIO.parse(DIR+gene+".clustal", "clustal"):
        for record in alignment:
            print(record.id)
            print(record.seq)
            if record.id in align_org.keys():
                align_org[record.id] += str(record.seq)
            else: align_org[record.id] = str(record.seq)
FASTA = ""
for org in organisms:
    FASTA += ">"+org + "\n" + align_org[org] + "\n"
DIR = "/media/alle/Bioinformatic/Bioinformatic/"
f = open(DIR+ "phylogenetictree" + ".fasta", "w")
f.write(FASTA)
f.close()

#creation of a table with each peptide hormone and motif score found for each organism used in the phylogenetic tree
DATABASE = psiBLAST_filtered3[1]
for PROPEPTIDE in DATABASE.keys():
    print(PROPEPTIDE)
    for ORGANISM in organisms:
        if ORGANISM in list(DATABASE[PROPEPTIDE]["ORG0"]):
            print("\t"+ORGANISM)

x= ""
for PROPEPTIDE in DATABASE.keys():
    x += "\t" + PROPEPTIDE
print(x)
for ORG in organisms:
    x = ORG
    for PROPEPTIDE in DATABASE.keys():
        if ORG in list(DATABASE[PROPEPTIDE]["ORG0"]):
            score = int(DATABASE[PROPEPTIDE][DATABASE[PROPEPTIDE]["ORG0"]==ORG]["motif_score"])
            x += "\t"+str(score)
        else: x += "\t"
    print(x)

        
# creation of the ensemble id to organism map
DIR = r"D:\Bioinformatic\blastdb\ensemble\zipped"
ORG0 = []
ORG1 = []
ORG2 = []
LIST = []
for filename in os.listdir(DIR):
    if filename.endswith(".fa"): 
        ORG = filename.split(".", 1)[0]
        print(os.path.join(DIR, filename))
        print(ORG)
        with open(os.path.join(DIR, filename)) as f:
            content = f.readlines()
            for x in content:
                if x.startswith(">"):
                    LIST.append(x.split(" ")[0][1:])
                    ORG0.append(ORG)
                    ORG1.append(ORG.split("_",1)[0])
                    ORG2.append(ORG.split("_",1)[1])

df = pd.DataFrame(data={
    'ACCESSION': LIST,
    'ORG0': ORG0,
    'ORG1': ORG1,
    'ORG2': ORG2
    })

df.sort_values(by=['ACCESSION'], inplace=True)
filename = str(DIR+"/Ensemble_accession_to_organism.csv")
df.to_csv(filename,
          sep = "\t",
          index=False
          )
                            
        