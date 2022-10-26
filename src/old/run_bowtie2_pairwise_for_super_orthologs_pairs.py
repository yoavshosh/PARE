# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 11:33:42 2019

@author: shosh
"""

import sys
import re
import os
import pandas as pd
import numpy as np
import itertools as it
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, IUPAC
import argparse
import time


def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['component','protein','orfs_start','orfs_end','strand','sequence']
    
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        if rec_data[0][-1] == ' ':  #some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(' ','')
        protein = rec_data[-1].split('|')[1]   #reading proteing from description assuming it was added to header using the transcriptome built pipeline we have for trinity
        rec_data = (rec_data[0],protein,int(rec_data[2]),int(rec_data[4]),rec_data[6],record.seq)
        data.append(rec_data)
    
    df = pd.DataFrame(data = data, columns = col_names)
    
    return df


def create_super_orthologs_df(all_best_hits,animals):
    """
    from all_best_hits dataframe containing rows of 2-way best hits orthologs for pairs of animals
    return a dataframe in which each row contain sequences for multiple animals that are all the 2-way best hits of eachother.
    """

    super_orthologs = None
    for pair in it.combinations(animals,2):
        a1 = pair[0]
        a2 = pair[1]
        pair_best_hits = all_best_hits[np.logical_or(np.logical_and(all_best_hits['animal_1']==a1,all_best_hits['animal_2']==a2),np.logical_and(all_best_hits['animal_1']==a2,all_best_hits['animal_2']==a1))]
        
        #swapping a2 and a2 (if needed) to fit correct order of columns
        a1 = pair_best_hits['animal_1'].iloc[0]
        a2 = pair_best_hits['animal_2'].iloc[0]

        #keeping only relevant columns in new names
        pair_best_hits = pair_best_hits[['score','seq1','seq2']]
        pair_best_hits = pair_best_hits.rename(index = str, columns = {'score':a1+'_'+a2+'_orthomcl_score','seq1':a1,'seq2':a2}) 
        
        if super_orthologs is None:
            super_orthologs = pair_best_hits
        else:
            if all(a in super_orthologs.columns for a in [a1,a2]):
#                print('both animals are already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = [a1 ,a2], how = 'inner')
            elif a1 in super_orthologs.columns and a2 not in super_orthologs.columns:
#                print('animal_1 (' + a1 + ') is already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = a1 , how = 'inner')
            elif a2 in super_orthologs.columns and a1 not in super_orthologs.columns:
#                print('animal_2 (' + a2 + ') is already in df')
                super_orthologs = super_orthologs.merge(pair_best_hits, on = a2, how = 'inner')
            if all(a not in super_orthologs.columns for a in [a1,a2]):
                print('NO INTERSECTION OF ANIMALS')
    
    return super_orthologs


if __name__ == '__main__':
    
    all_best_hits_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_8_species/all_best_hits.txt' 
    animals1 = ['oct','bim','sep','squ','bob','lin','nau']
    animals2 = ['oct','bim','sep','squ','bob','lin']
    animals3 = ['oct','bim','sep','squ','bob']
    animals4 = ['oct','bim','sep','squ','lin']
    animals5 = ['oct','bim','sep','squ']
    animals6 = ['apl','oct','bim','sep','squ','bob','lin','nau']
    animals7 = ['apl','oct','bim','sep','squ','nau']

    
    all_best_hits_file_name = all_best_hits_file.split('/')[-1]
    all_best_hits_file_path = '/'.join(all_best_hits_file.split('/')[:-1])+'/'
    
    print('Reading all best hits')
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    
    print('Creating super orthologs table')
    super_orthologs1 = create_super_orthologs_df(all_best_hits,animals1)
    super_orthologs2 = create_super_orthologs_df(all_best_hits,animals2)
    super_orthologs3 = create_super_orthologs_df(all_best_hits,animals3)
    super_orthologs4 = create_super_orthologs_df(all_best_hits,animals4)
    super_orthologs5 = create_super_orthologs_df(all_best_hits,animals5)
    super_orthologs6 = create_super_orthologs_df(all_best_hits,animals6)
    super_orthologs7 = create_super_orthologs_df(all_best_hits,animals7)
    
    
    