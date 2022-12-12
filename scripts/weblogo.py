# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 09:06:58 2022

@author: Daniel Moser AG Döhlemann Uni Cologne
"""
# Creating weblogo comparing AtPEP1 und ZmPEP1 motifs

# Import Packages
import pandas as pd
from Bio import AlignIO
import matplotlib.pyplot as plt
import logomaker
pd.options.mode.chained_assignment = None  # default='warn'

#File Path + Data path
PATH = "YOUR PATH/phytocytokine/align_AtPEP1_ZmPEP1/"
AtPEP1 = "AtPEP1_after_align_ZmPEP1 - Copy.clustal"
ZmPEP1 = "ZmPEP1_after_align_AtPEP1 - Copy.clustal"

# Function for reading in CLUSTAL
def read_clustal(FILEDIR):
    print(FILEDIR)
    align = AlignIO.read(FILEDIR, "clustal")
    print(align)
    LIST = []
    for record in align._records:
        LIST.append(str(record.seq))
    return(LIST)
       
# Function of creating Weblogo
def create_weblogo(FASTA, name):
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'bold',
            'size': 16,
            }

    df = {}
    df = logomaker.alignment_to_matrix(FASTA, to_type="counts")
    LOGO = logomaker.Logo(df,
                          color_scheme="dmslogo_funcgroup", #or chemistry
                          font_name="Times New Roman",
                          stack_order="big_on_top",
                          show_spines=True,
                          vpad=0.05)
    LOGO.style_xticks(rotation=0, fmt='%d', anchor=0, spacing=5)
    LOGO.draw(clear=False)
    plt.title(name, fontdict=font)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlabel("Position", fontdict=font)
    plt.ylabel("Count", fontdict=font)
    plt.tight_layout()
    plt.savefig(PATH+name+".png", dpi=600)
    plt.savefig(PATH+name+".pdf", dpi=600)
    plt.show()

# Execution
create_weblogo(read_clustal(PATH+AtPEP1), "AtPEP1")
create_weblogo(read_clustal(PATH+ZmPEP1), "ZmPEP1")
