# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 12:34:50 2019

@author: shosh

QA functions for orthomcl pipline
"""

import re
import os
import sys
import pandas as pd
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statsmodels.stats.multitest as p_adjust
from pylab import text
from scipy import stats
from collections import deque
from functools import reduce
from matplotlib import colors
from matplotlib.colors import LogNorm
from heapq import nsmallest
from Bio import SeqIO
import xlsxwriter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def check_intersection_of_pair_results_different_order(pair_df1,pair_df2,names = ['bob','lin']):    
    """
    for QA purposes - retuen all conserved sites shared by results from find_all_matching_mm_fixed.pl
    when executed with the parameters animal_1 animal_2
    and when executed with the parameters animal_2 animal_1
    """
    
    l1 = list(pair_df1.apply(lambda x: x[names[0]+'_key'] + '_' + x[names[1]+'_key'], axis = 1)) 
    l2 = list(pair_df2.apply(lambda x: x[names[0]+'_key'] + '_' + x[names[1]+'_key'], axis = 1))
    return list(set(l1)&set(l2))


def check_orthologs_scores(similar_sequences_path, seq1_id, seq2_id):
    
    columns = ['seq1','seq2','animal_1','animal_2','score1','score2','score3','score4']
    similar_sequences = pd.read_csv(similar_sequences_path, sep = '\t', names = columns)


def compare_sequences(query_fasta, subject_fasta):
    """
    check if sequences in query_fasta are all contained within subject_fasta by id
    and that they are all identical
    """
    comparison_path = '/'.join(query_fasta.split('/')[:-1])
    
    subject_fasta_dict = {}
    for srecord in SeqIO.parse(open(subject_fasta, "r"), "fasta"):
        subject_fasta_dict.update({srecord.id:str(srecord.seq)})
        
    query_records_not_in_subject = []
    query_records_not_identical_to_subject = []
    with open(comparison_path + '/comparison.txt', 'w') as comparison_file:
        for qrecord in SeqIO.parse(open(query_fasta, "r"), "fasta"):
            if qrecord.id in subject_fasta_dict:
                if str(qrecord.seq) != subject_fasta_dict[qrecord.id]:
                   query_records_not_identical_to_subject.append(qrecord)
                   comparison_file.write(qrecord.id + '\n')
                   comparison_file.write(str(qrecord.seq) + '\n')
                   comparison_file.write(subject_fasta_dict[qrecord.id] + '\n')  
            else:
                query_records_not_in_subject.append(qrecord.id)
    
    print(str(len(query_records_not_identical_to_subject)) + ' records not identical')
    print(str(len(query_records_not_in_subject)) + ' in ' + query_fasta + ' not in ' + subject_fasta)
    
    return query_records_not_in_subject, query_records_not_identical_to_subject

def find_lesser_hits_in_all_best_hits(all_best_hits,orthologs,similar_sequences):
    
    for index, row in all_best_hits.iterrows():
        seq = row['seq1']
        ortholog = row['seq2']
        ortholog_animal = ortholog[0:3]
        hits = orthologs[np.logical_or(orthologs['seq1']==seq,orthologs['seq2']==seq)]
        hits_for_pair = hits[np.logical_or(hits['animal_1']==ortholog_animal,hits['animal_2']==ortholog_animal)]
        
        similar_sequences_for_seq = similar_sequences[np.logical_or(similar_sequences['seq1']==seq,similar_sequences['seq2']==seq)]
        similar_sequences_for_seq = similar_sequences_for_seq[np.logical_or(similar_sequences_for_seq['animal_1']==ortholog_animal,similar_sequences_for_seq['animal_2']==ortholog_animal)]
        
        best_hit_1 = hits_for_pair.groupby('animal_1', group_keys=False).apply(lambda row: row.loc[row['score3'].idxmax()])
        best_hit_2 = hits_for_pair.groupby('animal_1', group_keys=False).apply(lambda row: row.loc[row['score4'].idxmax()])
        

def create_keys(row, k):

    name1 = k.split('_')[0]
    name2 = k.split('_')[1]
    try:    
        row['orthologs'] = row['Species #1 Trinity name'] + '_' + row['Species #2 Trinity name']
        row[name1+'_'+name2+'_'+'key'] = row['Species #1 Trinity name'] + '_'  + str(row['Species #1 Editing location inside sequence']) + '_' + row['Species #2 Trinity name'] + '_' + str(row['Species #2 Editing location inside sequence'])
    except KeyError:
        row['orthologs'] = row[name1+'_seq_id'] + '_' + row[name2+'_seq_id']
        row[name1+'_'+name2+'_'+'key'] = row[name1+'_seq_id'] + '_'  + str(row[name1+'_position']) + '_' + row[name2+'_seq_id'] + '_' + str(row[name2+'_position'])
    return row


def create_keys_all(row, animals):
    key = ''
    orthologs = ''
    try:
        for a in animals:
            key += '_' + row[a+' Trinity name'] + '_' + str(row[a+' Editing location inside sequence'])
            orthologs += '_' + row[a+' Trinity name']
    except KeyError:
        for a in animals:
            key += '_' + row[a+'_seq_id'] + '_' + str(row[a+'_position'])
            orthologs += '_' + row[a+'_seq_id']
    row['orthologs'] = orthologs[1:]
    row['key_all'] = key[1:]
    return row


def update_dfs_dict(dfs_dict, animals):
    updated_dict = {}
    for k, v in dfs_dict.items():
        if k != 'conserved_across_all':
            v = v.apply(lambda row: create_keys(row, k), axis = 1)
        else:
            v = v.apply(lambda row: create_keys_all(row, animals), axis = 1)
        updated_dict.update({k:v})
    return updated_dict


def compair_results(dfs_dict_1,dfs_dict_2, dfs_dict_2_name):
    
    dfs_dict_1_new = {}
    
    for k,v in dfs_dict_1.items():
        if k != 'conserved_across_all':
            keys_for_comparison = list(dfs_dict_2[k][k+'_key'])
            orthologs_for_comparison = list(dfs_dict_2[k]['orthologs'])
            v['key_in_'+dfs_dict_2_name] = v.apply(lambda row: row[k+'_key'] in keys_for_comparison, axis = 1)
            v['orthologs_in_'+dfs_dict_2_name] = v.apply(lambda row: row['orthologs'] in orthologs_for_comparison, axis = 1)        
        else:
            keys_for_comparison = list(dfs_dict_2[k]['key_all'])
            orthologs_for_comparison = list(dfs_dict_2[k]['orthologs'])
            v['key_in_'+dfs_dict_2_name] = v.apply(lambda row: row['key_all'] in keys_for_comparison, axis = 1)
            v['orthologs_in_'+dfs_dict_2_name] = v.apply(lambda row: row['orthologs'] in orthologs_for_comparison, axis = 1)        
        
        dfs_dict_1_new.update({k:v})
    
    return dfs_dict_1_new
        
        
def create_unmatching_keys_dict(dfs_dict, comparison_name):
    
    unmatching_keys_dict = {}
    unmatching_orthologs_dict = {}
    
    for k,v in dfs_dict.items():
        v_keys = v[~v['key_in_'+comparison_name]]
        unmatching_keys_dict.update({k:v_keys})
        print(k + ' keys not in ' + comparison_name + ': ' + str(len(v_keys)))        
        
        v_orthologs = v[~v['orthologs_in_'+comparison_name]]
        unmatching_orthologs_dict.update({k:v_orthologs})
        print(k + ' orthologs not in ' + comparison_name + ': ' + str(len(set(v_orthologs['orthologs']))))        
        
    return unmatching_keys_dict, unmatching_orthologs_dict


def all_best_hits_orthologs_keys(row, animals_order = ['bim','squ','oct','sep','bob','lin']):
    animal_1 = row['animal_1']
    animal_2 = row['animal_2']
    if animals_order.index(animal_1) < animals_order.index(animal_2):
        orthologs = row.seq1[4:] + '_' + row.seq2[4:]
    else:
        orthologs = row.seq2[4:] + '_' + row.seq1[4:]
    return orthologs
    



if __name__ == '__main__':
    
    animals = ['bim','oct','sep','squ']
    animals_noa = ['Oct.bim.','Oct.vul.','Sepia','Squid']

    all_best_hits_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/all_best_hits.txt'
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    all_best_hits['orthologs_key'] = all_best_hits.apply(lambda row: all_best_hits_orthologs_keys(row), axis=1)
    
    
    noa_old_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/results/analysis/conserved_old.xlsx'        
    noa_new_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/results/analysis/conserved_new.xlsx'        
    new_res_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/results/analysis/bim_squ_oct_sep.xlsx'    
    noa_old_res = pd.read_excel(noa_old_file, sheet_name=None)
    noa_new_res = pd.read_excel(noa_new_file, sheet_name=None)
    yoav_res = pd.read_excel(new_res_file, sheet_name=None)
    
    noa_old_res = update_dfs_dict(noa_old_res, animals_noa)
    noa_new_res = update_dfs_dict(noa_new_res, animals_noa)
    yoav_res = update_dfs_dict(yoav_res, animals)
    
    
    noa_old_res = compair_results(noa_old_res, noa_new_res, "noa_new")
    noa_old_res = compair_results(noa_old_res, yoav_res , "yoav")
    noa_new_res = compair_results(noa_new_res , noa_old_res , "noa_old")
    noa_new_res = compair_results(noa_new_res , yoav_res , "yoav")
    yoav_res = compair_results(yoav_res, noa_new_res, "noa_new")
    yoav_res = compair_results(yoav_res , noa_old_res , "noa_old")
    
    print('\nnoa_old - results not in other analyses')
    noa_old_not_in_noa_new = create_unmatching_keys_dict(noa_old_res, "noa_new")
    noa_old_not_in_yoav = create_unmatching_keys_dict(noa_old_res, "yoav")
    print('\nnoa_new - results not in other analyses')
    noa_new_not_in_noa_old = create_unmatching_keys_dict(noa_new_res, "noa_old")
    noa_new_not_in_yoav = create_unmatching_keys_dict(noa_new_res, "yoav")
    print('\nyoav - results not in other analyses')
    yoav_not_in_noa_old = create_unmatching_keys_dict(yoav_res, "noa_old")
    yoav_not_in_noa_new = create_unmatching_keys_dict(yoav_res, "noa_new")

    
    print('\norhologs in noa_new but not in last analysis best hits')
    for k,v in noa_new_not_in_yoav[1].items():
        orthologs_not_in_best_hits = set([i for i in list(v['orthologs']) if i not in list(all_best_hits['orthologs_key'])])
        print(k + ' ' + str(len(orthologs_not_in_best_hits)))
        print(orthologs_not_in_best_hits)
    
    print('\norhologs in noa_old but not in last analysis best hits')
    for k,v in noa_old_not_in_yoav[1].items():
        orthologs_not_in_best_hits = set([i for i in list(v['orthologs']) if i not in list(all_best_hits['orthologs_key'])])
        print(k + ' ' + str(len(orthologs_not_in_best_hits)))
        print(orthologs_not_in_best_hits)
    
    
    

    
###############################################################################


# =============================================================================
#     pair = 'oct_sep'
#     
#     name1 = pair.split('_')[0]
#     name2 = pair.split('_')[1]
#     
#     best_hits_for_pair = all_best_hits[np.logical_and(all_best_hits['animal_1']==name1, all_best_hits['animal_2']==name2)]
#     orthologs_from_best_hits = list(best_hits_for_pair['orthologs_key'])
#     
#     oct_sep_new = new_res_dfs[pair]
#     oct_sep_old = noa_old_res[pair]
#     new_oct_sep_not_in_old = new_res_not_in_old[pair]
#     old_oct_sep_not_in_new = noa_old_not_in_new[pair]
#     old_not_in_new_not_in_best_hits = [i for i in list(old_oct_sep_not_in_new['orthologs_key']) if i not in orthologs_from_best_hits]
#     new_not_in_old_not_in_best_hits = [i for i in list(new_oct_sep_not_in_old['orthologs_key']) if i not in orthologs_from_best_hits]
#     
#     #compairing the orthologs that were not mapped correctly to orthologs in the other list (not including the specific sites mapping in potential orthologs)
#     #this is just to verify that wrong pairing is in proteins dataset level and not just different mapping of sites in a correctly maped protein 
#     new_was_not_mapped_as_old = [i for i in list(new_oct_sep_not_in_old['orthologs_key']) if i not in list(noa_old_res[pair]['orthologs_key'])]
#     old_was_not_mapped_as_new = [i for i in list(old_oct_sep_not_in_new['orthologs_key']) if i not in list(new_res_dfs[pair]['orthologs_key'])]
#     
#     new_mapped_as_old = [i for i in list(new_oct_sep_not_in_old['orthologs_key']) if i in list(noa_old_res[pair]['orthologs_key'])]
#     old_mapped_as_new = [i for i in list(old_oct_sep_not_in_new['orthologs_key']) if i in list(new_res_dfs[pair]['orthologs_key'])]
# =============================================================================
    
    
###############################################################################
    #reading all potential orthologs and all similarSequences data from blast stages of orthomcl
# =============================================================================
#     similarSequences_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/similarSequences.txt'
#     similarSequences_file = pd.read_csv(similarSequences_file, sep = '\t', names = ['seq1','seq2','animal_1','animal_2','evalueMant','evalueExp','identity','match'])
#     
#     orthologs_file = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/orthologs.txt'
#     orthologs = pd.read_csv(orthologs_file, sep = '\t', names = ['seq1','seq2','score'])  
#     orthologs['animal_1'] = orthologs.apply(lambda row: row['seq1'][0:3], axis = 1)
#     orthologs['animal_2'] = orthologs.apply(lambda row: row['seq2'][0:3], axis = 1)
# =============================================================================


"""
noa_old - results not in other analyses
bim_oct keys not in noa_new: 73
bim_oct orthologs not in noa_new: 7
squ_sep keys not in noa_new: 75
squ_sep orthologs not in noa_new: 7
bim_sep keys not in noa_new: 30
bim_sep orthologs not in noa_new: 10
squ_oct keys not in noa_new: 40
squ_oct orthologs not in noa_new: 5
oct_sep keys not in noa_new: 69
oct_sep orthologs not in noa_new: 14
bim_squ keys not in noa_new: 16
bim_squ orthologs not in noa_new: 4
conserved_across_all keys not in noa_new: 5
conserved_across_all orthologs not in noa_new: 1
bim_oct keys not in yoav: 408
bim_oct orthologs not in yoav: 67
squ_sep keys not in yoav: 473
squ_sep orthologs not in yoav: 52
bim_sep keys not in yoav: 135
bim_sep orthologs not in yoav: 27
squ_oct keys not in yoav: 125
squ_oct orthologs not in yoav: 22
oct_sep keys not in yoav: 160
oct_sep orthologs not in yoav: 38
bim_squ keys not in yoav: 95
bim_squ orthologs not in yoav: 26
conserved_across_all keys not in yoav: 109
conserved_across_all orthologs not in yoav: 48

noa_new - results not in other analyses
bim_oct keys not in noa_old: 60
bim_oct orthologs not in noa_old: 0
squ_sep keys not in noa_old: 204
squ_sep orthologs not in noa_old: 12
bim_sep keys not in noa_old: 183
bim_sep orthologs not in noa_old: 33
squ_oct keys not in noa_old: 192
squ_oct orthologs not in noa_old: 26
oct_sep keys not in noa_old: 292
oct_sep orthologs not in noa_old: 46
bim_squ keys not in noa_old: 139
bim_squ orthologs not in noa_old: 27
conserved_across_all keys not in noa_old: 64
conserved_across_all orthologs not in noa_old: 16
bim_oct keys not in yoav: 337
bim_oct orthologs not in yoav: 60
squ_sep keys not in yoav: 402
squ_sep orthologs not in yoav: 47
bim_sep keys not in yoav: 109
bim_sep orthologs not in yoav: 20
squ_oct keys not in yoav: 87
squ_oct orthologs not in yoav: 17
oct_sep keys not in yoav: 94
oct_sep orthologs not in yoav: 26
bim_squ keys not in yoav: 84
bim_squ orthologs not in yoav: 23
conserved_across_all keys not in yoav: 113
conserved_across_all orthologs not in yoav: 50

yoav - results not in other analyses
bim_squ keys not in noa_old: 519
bim_squ orthologs not in noa_old: 90
bim_oct keys not in noa_old: 538
bim_oct orthologs not in noa_old: 7
bim_sep keys not in noa_old: 662
bim_sep orthologs not in noa_old: 111
squ_oct keys not in noa_old: 706
squ_oct orthologs not in noa_old: 107
squ_sep keys not in noa_old: 874
squ_sep orthologs not in noa_old: 49
oct_sep keys not in noa_old: 921
oct_sep orthologs not in noa_old: 124
conserved_across_all keys not in noa_old: 202
conserved_across_all orthologs not in noa_old: 32
bim_squ keys not in noa_new: 385
bim_squ orthologs not in noa_new: 64
bim_oct keys not in noa_new: 480
bim_oct orthologs not in noa_new: 7
bim_sep keys not in noa_new: 483
bim_sep orthologs not in noa_new: 81
squ_oct keys not in noa_new: 516
squ_oct orthologs not in noa_new: 81
squ_sep keys not in noa_new: 674
squ_sep orthologs not in noa_new: 39
oct_sep keys not in noa_new: 632
oct_sep orthologs not in noa_new: 80
conserved_across_all keys not in noa_new: 147
conserved_across_all orthologs not in noa_new: 19

orhologs in noa_new but not in last analysis best hits
bim_oct 8
squ_sep 9
bim_sep 3
squ_oct 2
oct_sep 2
bim_squ 5
conserved_across_all 50

orhologs in noa_old but not in last analysis best hits
bim_oct 8
squ_sep 9
bim_sep 3
squ_oct 2
oct_sep 2
bim_squ 5
conserved_across_all 48
"""



