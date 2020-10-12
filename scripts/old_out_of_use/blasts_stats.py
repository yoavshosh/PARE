# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:38:27 2019

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


def read_blast_results(blast_file):
    """
    read blast resutls file (all-vs_all) and keep only records for a1 as query animal and a2 as subject animal
    """    
    data = []
#    col_names = ['comp_query', 'comp_subject','identity','alignment_length','mismatch','gapopen',
#               'query_start','query_end','subject_start','subject_end','e_value','bitscore','btop']
#    
    col_names = ['comp_query', 'comp_subject','identity','alignment_length','mismatch','gapopen',
               'query_start','query_end','subject_start','subject_end','e_value','bitscore','qcovhsp','qcovus','btop']

    
    with open(blast_file, "r") as bf:
        for i, line in enumerate(bf.readlines()):
            fields = line.split('\t')
            fields[-1] = fields[-1].replace('\n','')
            data.append(fields)
    blast_df = pd.DataFrame(data = data, columns = col_names)
    
    for c in col_names:
        blast_df[c] = pd.to_numeric(blast_df[c],errors = 'ignore')
    
    return blast_df

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


def grouped_bar_chart(pairs_stats_df, fig_path, values_column, labels_column = 'pair'):
    
    index = np.arange(len(pairs_stats_df))
    plt.bar(index, pairs_stats_df[values_column])
    plt.title('Super Orthologs Protein Pairs Average Bitscore', fontsize=10)
    plt.xlabel('Pair', fontsize=10)
    plt.ylabel('Average bitscore', fontsize=10)
    plt.xticks(index, pairs_stats_df[labels_column], fontsize=10, rotation=90)
    plt.ylim(0,max(pairs_stats_df[values_column]*1.1))
    
    plt.savefig(fig_path+values_column+'.jpg',bbox_inches='tight')
#    plt.margins(x=1.3, y=1.3)
#    plt.show()
    plt.close()
    
 
def collect_super_orthologs_blast(blastp, super_orthologs, animals):
    """
    This function creates a blasts results file contaeining only results of pairs that are super orthologs
    """
    top_scores_blasts = []
    for a1,a2 in it.combinations(animals, 2):
        
        animals_blast = blastp[np.logical_and(blastp['query_animal'].isin([a1,a2]),blastp['subject_animal'].isin([a1,a2]))]
        
        for seq1,seq2 in zip(super_orthologs[a1],super_orthologs[a2]):
            pair = [seq1,seq2]
            sequences_blast = animals_blast[np.logical_or(np.logical_and(animals_blast['comp_query']==seq1,animals_blast['comp_subject']==seq2),np.logical_and(animals_blast['comp_query']==seq2,animals_blast['comp_subject']==seq1))]
            try:
                sequences_blast = sequences_blast.loc[sequences_blast['bitscore'].idxmax(),:]
                top_scores_blasts.append(sequences_blast)
            except ValueError:
                pass
            
    super_orthologs_blast = pd.DataFrame(top_scores_blasts)
    
    return super_orthologs_blast        
          

def plot_bitscores_std_hist(bitscores_std, all_best_hits_file_path, bins='auto'):
    
    n, bins, patches = plt.hist(x=bitscores_std, bins=bins, color='#0504aa', alpha=0.7, rwidth=0.85, density=True)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('super - orthologs bitscore std')
    plt.ylabel('Frequency')
#    plt.text(23, 45, r'$\mu=15, b=3$')
#    maxfreq = n.max()
# Set a clean upper y-axis limit.
#    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(all_best_hits_file_path+'super_orthologs_bitscores_std.jpg',bbox_inches='tight')
    plt.close()


def collect_super_orthologs_bit_scores(row,super_orthologs_blast,animals):
    """
    Collect bitscores from blast results for each pair in super orthologs group
    """    
    for a1,a2 in it.combinations(animals, 2):
        row[a1+'_'+a2+'_bitscore'] = super_orthologs_blast[np.logical_or(np.logical_and(super_orthologs_blast['comp_query']==row[a1],super_orthologs_blast['comp_subject']==row[a2]),np.logical_and(super_orthologs_blast['comp_query']==row[a2],super_orthologs_blast['comp_subject']==row[a1]))].squeeze()['bitscore']
    return row


def calculate_pairs_distribution_per_gene(row):
    average_bitscore = np.mean(row)    
    for p in row.index:
        row[p+'_distance'] = row[p] - average_bitscore
    row['average_bitscore'] = average_bitscore
    row['bitscores_std'] = np.std(row.values)
    return row
    

def plot_distances_distributions(bitscore_df, path):
    
    palette = plt.get_cmap('Set1')
    columns = [c for c in bitscore_df.columns if not(any([s in c for s in ('_distance','_std')]))]
    for i,col in enumerate(columns):
        if i==0:
            palette = plt.get_cmap('Set1')
        elif 'nau' in col:
            palette = plt.get_cmap('plasma')
        else:
            palette = plt.get_cmap('spring_r')
        plt.plot(list(bitscore_df.index), bitscore_df[col], marker='', color=palette(i), linewidth=1, alpha=1, label=col)

    plt.legend(loc=(1,0) , ncol=2)
    plt.xlabel("gene index")
    plt.ylabel("Pair difference from average")
    plt.savefig(all_best_hits_file_path+'distances_from_average_spaghetti.jpg',bbox_inches='tight')
    plt.close()


def all_pairs_blast(row, animals):
    """
    check if all pairs in row have blast results
    """
    all_pairs_have_blast = True
    for a1,a2 in it.combinations(animals, 2):
        pair = a1+'_'+a2
#        print(str(type(row[pair+'_bitscore'])))
        if type(row[pair+'_bitscore'])!=float:
            all_pairs_have_blast = False
            break
     
    return all_pairs_have_blast
          
    
if __name__ == '__main__':
    
    all_best_hits_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_8_species/blast_stats/all_best_hits.txt' 
    blast_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_8_species/blast_stats/super_orthologs_blast.txt' 
    create_plot = True
    
#    all_best_hits_file = sys.argv[1] 
#    blast_file = sys.argv[2]
#    create_plot = eval(sys.argv[3])
    
    all_best_hits_file_name = all_best_hits_file.split('/')[-1]
    all_best_hits_file_path = '/'.join(all_best_hits_file.split('/')[:-1])+'/'
    
#    sys.stdout = open(all_best_hits_file_path+'blast_stats_super_orthologs_from_'+all_best_hits_file_name, 'w')
    
    animals = ['oct','bim','sep','squ','bob','lin','nau']
    
    print('Reading all best hits')
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    
    print('Creating super orthologs table')
    super_orthologs = create_super_orthologs_df(all_best_hits,animals)
    
    print('Reading Blast results')
    blastp = read_blast_results(blast_file)
    blastp['query_animal'] = blastp.apply(lambda row: row['comp_query'].split('|')[0], axis = 1)
    blastp['subject_animal'] = blastp.apply(lambda row: row['comp_subject'].split('|')[0], axis = 1)
    
    print('Collecting super orthologs blast results')
    super_orthologs_blast = collect_super_orthologs_blast(blastp, super_orthologs, animals)
    super_orthologs_blast.to_csv(all_best_hits_file_path+'super_orthologs_blast.txt', sep='\t', index=False)
    super_orthologs_blast = pd.read_csv(all_best_hits_file_path+'super_orthologs_blast.txt', sep='\t')
    
    print('Collecting super orthologs bistscores to super orthologs table')
    super_orthologs = super_orthologs.apply(lambda row: collect_super_orthologs_bit_scores(row,super_orthologs_blast,animals), axis=1)           
    super_orthologs['all_have_blast'] =  super_orthologs.apply(lambda row: all_pairs_blast(row, animals), axis=1) 
    super_orthologs = super_orthologs[super_orthologs['all_have_blast']]
    super_orthologs.to_csv(all_best_hits_file_path+'super_orthologs_scores.txt', sep='\t', index=False)

    print('Pairs stats:')
    columns = ['pair','orthomcl_score','blastp_bitscore']
    stats = []
    for a1,a2 in it.combinations(animals, 2):
        pair = a1+'_'+a2
        print(pair)
        try:
            orthomcl_score_mean = super_orthologs[a1+'_'+a2+'_orthomcl_score'].mean()
        except KeyError:
            orthomcl_score_mean = super_orthologs[a2+'_'+a1+'_orthomcl_score'].mean()
        print('Average orthomcl super orthologs score (from all_best_hits): ' + str(orthomcl_score_mean))    
        
        try:
            bitscore_mean = super_orthologs[a1+'_'+a2+'_bitscore'].mean()
        except KeyError:
            bitscore_mean = super_orthologs[a2+'_'+a1+'_bitscore'].mean()
        print('Average super orthologs bitscore: ' + str(bitscore_mean))
        stats.append([pair,orthomcl_score_mean,bitscore_mean])
    pairs_stats_df = pd.DataFrame(columns=columns,data=stats)
    
    super_orthologs_blast = pd.read_csv(all_best_hits_file_path+'super_orthologs_blast.txt', sep='\t')
    super_orthologs = pd.read_csv(all_best_hits_file_path+'super_orthologs_scores.txt', sep='\t')
    if create_plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        bitscore_columns = [c for c in list(super_orthologs.columns) if 'bitscore' in c]
        bitscore_df = super_orthologs.loc[:,bitscore_columns]
        
        for col in list(bitscore_df.columns):
            bitscore_df = bitscore_df.rename(columns={col: col.replace('_bitscore','')})  
        data=pd.melt(bitscore_df)
        data['value']=data['value'].astype(float)
        data = data.rename(columns={"value": "Bitscore", "variable": "Pair"})
        g = sns.boxplot(x="Pair", y="Bitscore", data=data, fliersize = 0.5)
#        g = sns.violinplot(x="Pair", y="Bitscore", data=data, fliersize = 0.5)
        plt.xticks(rotation=90)
        plt.savefig(all_best_hits_file_path+'bitscores_boxplot.jpg',bbox_inches='tight')
        plt.close()
            
        bitscore_df = bitscore_df.apply(lambda row: calculate_pairs_distribution_per_gene(row), axis=1)
        np.histogram(bitscore_df['bitscores_std'], bins=100)
        
        bitscore_dist_columns = [c for c in list(bitscore_df.columns) if '_distance' in c]
        bitscore_dist_df = bitscore_df.loc[:,bitscore_dist_columns]
        for col in list(bitscore_dist_df.columns):
            bitscore_dist_df = bitscore_dist_df.rename(columns={col: col.replace('_distance','')})  
      
        plot_distances_distributions(bitscore_dist_df,all_best_hits_file_path)
          
        data=pd.melt(bitscore_dist_df)
        data['value']=data['value'].astype(float)
        data ['Pair'] = data.apply(lambda row: row["variable"].replace('_','-'), axis=1)
        data = data.rename(columns={"value": "Distances from average bitscore"})
        g = sns.boxplot(x="Pair", y="Distances from average bitscore", data=data, fliersize = 0.5)
        
#        g = sns.violinplot(x="Pair", y="Distances from average bitscore", data=data, fliersize = 0.5)
        plt.xticks(rotation=90)
        plt.savefig(all_best_hits_file_path+'distances_from_average_bitscore_boxplot.jpg',bbox_inches='tight')
        plt.close()
        
        grouped_bar_chart(pairs_stats_df, all_best_hits_file_path, 'orthomcl_score', labels_column = 'pair')
        grouped_bar_chart(pairs_stats_df, all_best_hits_file_path, 'blastp_bitscore', labels_column = 'pair')
        
   
# =============================================================================
#     print('Pairs stats:')
#     columns = ['pair','orthomcl_score','blastp_bitscore']
#     stats = []
#     for a1,a2 in it.combinations(animals, 2):
#         pair = a1+'_'+a2
#         print(pair)
#         animals_blast = blastp[np.logical_and(blastp['query_animal'].isin([a1,a2]),blastp['subject_animal'].isin([a1,a2]))]
#         try:
#             score_mean = super_orthologs[a1+'_'+a2+'_score'].mean()
#         except KeyError:
#             score_mean = super_orthologs[a2+'_'+a1+'_score'].mean()
#         print('Average orthomcl super orthologs score (from all_best_hits): ' + str(score_mean))
#         
#         bitscore_sum=0
#         for seq1,seq2 in zip(super_orthologs[a1],super_orthologs[a2]):
#             pair = [seq1,seq2]
#             sequences_blast = animals_blast[np.logical_and(blastp['comp_query'].isin(pair),blastp['comp_subject'].isin(pair))]
#             if len(sequences_blast) > 2:
#                 print(pair)
#                 print('more than 2 blast results for pair')
#             sequences_blast = sequences_blast.loc[sequences_blast['bitscore'].idxmax(),:]        
#             bitscore_sum+=sequences_blast['bitscore']
#             average_bitscore = bitscore_sum/len(super_orthologs)
#         print('Average super orthologs bitscore: ' + str(average_bitscore))
#         stats.append([pair,score_mean,average_bitscore])
#     
#     pairs_stats_df = pd.DataFrame(columns=columns,data=stats)
# 
#     if create_plot:
#         import matplotlib.pyplot as plt
#         fig_path = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/'
#         grouped_bar_chart(pairs_stats_df, fig_path, 'orthomcl_score', labels_column = 'pair')
#         grouped_bar_chart(pairs_stats_df, fig_path, 'blastp_bitscore', labels_column = 'pair')
# =============================================================================

# =============================================================================
#     #this block is good for reconstructing the stats df based on stdout
#     columns = ['pair','orthomcl_score','blastp_bitscore']
#     stats = []
#     pair_data = []
#     for i,line in enumerate(text.split('\r\n')):
#         if i%3==0:
#             pair_data = []    
#             pair_data.append(line)
#         elif i%3==1:
#             pair_data.append(float(line.split(' ')[-1]))
#         elif i%3==2:
#             pair_data.append(float(line.split(' ')[-1]))
#             stats.append(pair_data)
# =============================================================================
   

"""
oct_bim
Average super orthologs score from all_best_hits: 1.0566614567052206
Average super orthologs bitscore: 1210.2585182869648
oct_sep
Average super orthologs score from all_best_hits: 1.0787977492966554
Average super orthologs bitscore: 1232.747733666771
oct_squ
Average super orthologs score from all_best_hits: 1.0459365426695844
Average super orthologs bitscore: 1239.6658330728353
oct_bob
Average super orthologs score from all_best_hits: 1.038740856517662
Average super orthologs bitscore: 1246.6195686151923
oct_lin
Average super orthologs score from all_best_hits: 1.043709909346671
Average super orthologs bitscore: 1239.423569865583
oct_nau
Average super orthologs score from all_best_hits: 1.0446398874648328
Average super orthologs bitscore: 1233.2100656455143
bim_sep
Average super orthologs score from all_best_hits: 1.0567308533916848
Average super orthologs bitscore: 1238.0678336980307
bim_squ
Average super orthologs score from all_best_hits: 1.0485282900906532
Average super orthologs bitscore: 1243.8859018443263
bim_bob
Average super orthologs score from all_best_hits: 1.0437086589559237
Average super orthologs bitscore: 1250.9571741169116
bim_lin
Average super orthologs score from all_best_hits: 1.0433104095029697
Average super orthologs bitscore: 1243.1981869334168
bim_nau
Average super orthologs score from all_best_hits: 1.023111909971866
Average super orthologs bitscore: 1236.9134104407628
sep_squ
Average super orthologs score from all_best_hits: 1.0533107221006563
Average super orthologs bitscore: 1230.4101281650517
sep_bob
Average super orthologs score from all_best_hits: 1.0415514223194748
Average super orthologs bitscore: 1239.6820881525477
sep_lin
Average super orthologs score from all_best_hits: 1.0404623319787432
Average super orthologs bitscore: 1230.5014066895906
sep_nau
Average super orthologs score from all_best_hits: 1.0601309784307595
Average super orthologs bitscore: 1229.2585182869648
squ_bob
Average super orthologs score from all_best_hits: 1.0247561738043136
Average super orthologs bitscore: 1244.137230384495
squ_lin
Average super orthologs score from all_best_hits: 1.0251225382932165
Average super orthologs bitscore: 1239.5154735854956
squ_nau
Average super orthologs score from all_best_hits: 1.0508446389496717
Average super orthologs bitscore: 1241.8208815254768
bob_lin
Average super orthologs score from all_best_hits: 1.0167339793685526
Average super orthologs bitscore: 1243.353235386058
bob_nau
Average super orthologs score from all_best_hits: 1.0528587058455767
Average super orthologs bitscore: 1248.3851203501094
lin_nau
Average super orthologs score from all_best_hits: 1.0406595811190997
Average super orthologs bitscore: 1234.9756173804315
"""
        
    
    
    