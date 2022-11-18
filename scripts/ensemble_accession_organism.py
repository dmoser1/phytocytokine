# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 18:45:07 2020

@author: dan2m
"""

# Author: Daniel Moser, AG Doehlemann, University of Cologne, Germany
# Purpose: Linking protein IDs to source organisms from plant proteoms retrieved from ensemble plants

# Loading packages
import os
import pandas as pd

# File directory
DIR = r"ENTER YOUR FILEPATH TO THE PROTEOMS (FASTA) HERE"

# Matching each ID with the organism name which was in the file name
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
                    