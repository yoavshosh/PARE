# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:56:27 2019

@author: shosh
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
import seaborn as sns
from pylab import text
from scipy import stats
from collections import deque
from functools import reduce
from matplotlib import colors
from matplotlib.colors import LogNorm
from heapq import nsmallest
from Bio import SeqIO
import xlsxwriter


def scatter_orthologs_editing_levels(output_folder,el_1,el_2,animal_1,animal_2):

    pearson_corel_coef, p_val = stats.pearsonr(el_1, el_2)
    fig, ax = plt.subplots()
    plt.plot(el_1,el_2, c = 'red', linestyle='', marker = 'o', markersize = '0.7')
    plt.xlabel(animal_1)
    plt.ylabel(animal_2)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks([0,0.5,1])
    plt.yticks([0,0.5,1])
    plt.autoscale()
    plt.title('Ortholog Editing sites -\nEditing Levels Comparison')
    plt.text(0.05,1,'Rho = '+str(round(pearson_corel_coef,2)))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
#    plt.show()
    plt.savefig(output_folder+animal_1+'_'+animal_2+'_editing_levels_of_conserved_sites.png')
    plt.close()


def fmt(x, pos):
    return r"10^{}".format(x)


if __name__ == '__main__':
   
    results_path = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/results/'
    file = 'all_sites_mappings.txt'
    animals = ['oct','bim','sep','squ','lin','bob','nau']
   
#    file = 'all_sites_mappings4.txt'
#    animals = ['oct','bim','sep','squ']
   
    if not os.path.exists(results_path + 'graphs/'+'_'.join(animals)+'/'):
        os.makedirs(results_path + 'graphs/'+'_'.join(animals)+'/')
    output_folder = results_path + 'graphs/'+'_'.join(animals)+'/'    

    columns_suffix = ('component','protein','edited','coding_location','nuc','loc_inside_codon','aa_change','site_key','editing_level')
    conserved_pairs = []
    conserved_sites = []
    conserved_recoding_sites = []
#    pairs_dict = {}

    grand_locations_table = pd.read_csv(results_path+file, sep='\t', index_col = False)
    excel_writer = pd.ExcelWriter(results_path+'_'.join(animals)+'.xlsx',engine='xlsxwriter')
    grand_locations_table.to_excel(excel_writer, sheet_name='all_mappings', startrow=0 , startcol=0, index=False) 
    for pair in it.combinations(animals,2):
        a1 = pair[0]
        a2 = pair[1]
        pair_columns = [pair[i]+'_'+col for i in range(len(pair)) for col in columns_suffix]
#        print(pair)
#        print(pair_columns)
        pair_df = grand_locations_table[pair_columns].copy()
#        pair_df = pair_df[np.logical_or(pair_df[a1+'_edited'],pair_df[a2+'_edited'])]
#        pairs_dict.update({a1+'_'+a2:pair_df})
        pair_df = pair_df[np.logical_and(pair_df[a1+'_edited'],pair_df[a2+'_edited'])]
        pair_df[a1+'_'+a2+'_conserved'] = pair_df.apply(lambda row: row[a1+'_aa_change'][1] == row[a2+'_aa_change'][1], axis = 1) #defining conserved sites as identical target AA
        pair_df.to_excel(excel_writer, sheet_name=a1+'_'+a2+'_es', startrow=0 , startcol=0, index=False)
       
        pair_df = pair_df[pair_df[a1+'_'+a2+'_conserved']]
        conserved_sites.append(len(pair_df))
        conserved_pairs.append(a1+'_'+a2)

        pair_df[a1+'_recoding'] = pair_df.apply(lambda row: row[a1+'_aa_change'][0]!=row[a1+'_aa_change'][1], axis = 1)
        pair_df[a2+'_recoding'] = pair_df.apply(lambda row: row[a2+'_aa_change'][0]!=row[a2+'_aa_change'][1], axis = 1)
       
    
        #keep only conserved edited sites that are recoding in both animals
        pair_df = pair_df[np.logical_and(pair_df[a1+'_recoding'],pair_df[a2+'_recoding'])]
        conserved_recoding_sites.append(len(pair_df))
       
        g = (sns.jointplot(pair_df[a1+'_editing_level'], pair_df[a2+'_editing_level'], kind="hex",bins='log',joint_kws = {'gridsize':30}).set_axis_labels(a1, a2))
#        g = (sns.jointplot(conserved_across_all[a1+'_editing_level'], conserved_across_all[a2+'_editing_level'], kind="hex",bins='log',joint_kws = {'gridsize':30}).set_axis_labels(a1, a2))
           
        cbar_ax = g.fig.add_axes([1, .25, .05, .4])
        plt.colorbar(cax=cbar_ax, format=ticker.FuncFormatter(fmt))
        g.savefig(output_folder+a1+'_'+a2+'_recoding_editing_levels_for_conserved_sites')
        plt.close()

    excel_writer.save()
    
    plt.bar(np.arange(len(conserved_sites)), conserved_sites, align='center', alpha=0.5)
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.xticks(np.arange(len(conserved_sites)), conserved_pairs, rotation=70)
    plt.ylabel('Conserved ES')
    plt.savefig(output_folder+'Conserved_sites.png') 
    plt.close()

    plt.bar(np.arange(len(conserved_recoding_sites)), conserved_recoding_sites, align='center', alpha=0.5)
    plt.gcf().subplots_adjust(bottom=0.2)
    plt.xticks(np.arange(len(conserved_recoding_sites)), conserved_pairs, rotation=70)
    plt.ylabel('Conserved Recoding events')
    plt.savefig(output_folder+'Conserved_Recoding_events.png') 
    plt.close()
   
    grand_locations_table['editing_events'] = grand_locations_table.apply(lambda row: sum([row[a+'_edited'] for a in animals]), axis = 1)
    for i in set(list(grand_locations_table['editing_events'])):
        #initialize animals dict
        
        editing_events_i = grand_locations_table[grand_locations_table['editing_events']==i]
        a_freq = []
        g_freq = []
        c_freq = []
        t_freq = []
        gap_freq = []
        nuc_labels = ['A','G','C','T','gap']

        unedited_locations_number = []
        animals_in_plot = []
        for a in animals:
           
            unedited_locations = editing_events_i[~editing_events_i[a+'_edited']]
            
            if len(unedited_locations):
                
                animals_in_plot.append(a)
                unedited_locations_number.append(len(unedited_locations))
                a_freq.append(len(unedited_locations[unedited_locations[a+'_nuc']=='A'])/len(unedited_locations))
                g_freq.append(len(unedited_locations[unedited_locations[a+'_nuc']=='G'])/len(unedited_locations))
                c_freq.append(len(unedited_locations[unedited_locations[a+'_nuc']=='C'])/len(unedited_locations))
                t_freq.append(len(unedited_locations[unedited_locations[a+'_nuc']=='T'])/len(unedited_locations))
                gap_freq.append(len(unedited_locations[unedited_locations[a+'_nuc']=='gap'])/len(unedited_locations))
                
        # Set position of bar on X axis
        barWidth = 0.15
        r1 = np.arange(len(animals_in_plot))
        r2 = [x + barWidth for x in r1]
        r3 = [x + barWidth for x in r2]
        r4 = [x + barWidth for x in r3]
        r5 = [x + barWidth for x in r4]
    
    
        # Make the plot
        plt.bar(r1, a_freq, color='b', width=barWidth, edgecolor='white', label='A')
        plt.bar(r2, g_freq, color='y', width=barWidth, edgecolor='white', label='G')
        plt.bar(r3, c_freq, color='g', width=barWidth, edgecolor='white', label='C')
        plt.bar(r4, t_freq, color='r', width=barWidth, edgecolor='white', label='T')
        plt.bar(r5, gap_freq, color='black', width=barWidth, edgecolor='white', label='gap')
        
        # Add xticks on the middle of the group bars
        plt.gcf().subplots_adjust(bottom=0.5)
        plt.xlabel('animals', fontweight='bold')
        plt.xticks([r + barWidth for r in range(len(animals_in_plot))], [a+' n='+str(unedited_locations_number[i]) for i,a in enumerate(animals_in_plot)],rotation=-90)
        plt.ylabel('Unedited nucleotides frequency')
        plt.title('Unedited Nucleotide Frequencies Where Edited Animals = ' + str(i))
        # Create legend & Show graphic
        plt.legend(bbox_to_anchor=(1.1, 1))
        plt.savefig(output_folder + 'unedited_freq_edited='+str(i)+'.png')
        plt.close()
        
# =============================================================================
# 
# def conserved_and_recoding_across_all(row, animals):
#     
#     target_aa = []
#     recoding = []
#     for a in animals:
#         target_aa.append(row[a+'_aa_change'][1])
#         if row[a+'_aa_change'][0]!=row[a+'_aa_change'][1]:
#             recoding.append(True)
#         else:
#             recoding.append(False)
#     row['recoding_across_all'] = all(recoding)
#     
#     if len(set(target_aa)) == 1:
#         row['conserved_across_all'] = True
#     else:
#         row['conserved_across_all'] = False
#     
#     return row
#     
# animals = ['oct','bim','sep','squ','lin','bob']
# 
# edited_across_all = grand_locations_table
# for a in animals:
#     edited_across_all = edited_across_all[edited_across_all[a+'_edited']]
#  
# edited_across_all = edited_across_all.apply(lambda row: conserved_and_recoding_across_all(row, animals), axis = 1)
# 
# conserved_across_all = edited_across_all[edited_across_all['conserved_across_all']]
# conserved_recoding_across_all = conserved_across_all[conserved_across_all['recoding_across_all']].copy()
# 
# for a in animals:
#     es_df = super_ortholog_editing_sites_dict[a]
#     average_coverage = round(np.mean(es_df['RNA_coverage']),2)
#     average_highly_conserved_coverage = round(np.mean(es_df[es_df['component'].isin(list(conserved_across_all[a+'_component']))]['RNA_coverage']),2)
#     print(a + ' coverage for sites in super orthologs: '+ str(average_coverage))
#     print(a + ' coverage for highly conserved sites: '+ str(average_highly_conserved_coverage))
#  
# conserved_recoding_across_all['average_editing_for_5'] =  conserved_recoding_across_all.apply(lambda row: sum([row['oct_editing_level'],row['bim_editing_level'],row['sep_editing_level'],row['squ_editing_level'],row['lin_editing_level']])/5, axis = 1)
# conserved_recoding_across_all['average_editing'] =  conserved_recoding_across_all.apply(lambda row: sum([row['oct_editing_level'],row['bim_editing_level'],row['sep_editing_level'],row['squ_editing_level'],row['lin_editing_level'],row['bob_editing_level']])/6, axis = 1)
# 
# conserved_recoding_across_all.sort_values(by='average_editing_for_5', inplace = True)
# x = np.arange(len(conserved_recoding_across_all))
#     
# plt.gca().set_color_cycle(['orange', 'blue','yellow','green','blue','black'])
# plt.plot(x,conserved_recoding_across_all['oct_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=0.3)
# plt.plot(x,conserved_recoding_across_all['bim_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=0.3)
# plt.plot(x,conserved_recoding_across_all['sep_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=0.3)
# plt.plot(x,conserved_recoding_across_all['squ_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=0.3)
# plt.plot(x,conserved_recoding_across_all['lin_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=1)
# #    plt.plot(x,conserved_recoding_across_all['average_editing_for_5'],linewidth=0.1)
# plt.plot(x,conserved_recoding_across_all['bob_editing_level']-conserved_recoding_across_all['average_editing'],linewidth=1)
# plt.xlabel('site index')
# plt.ylabel('deviation from average editin level')
# #    plt.legend(['oct_editing_level','bim_editing_level','sep_editing_level','squ_editing_level','lin_editing_level','average for 5 animals', 'bob editing level'], loc='upper left')
# plt.legend(['oct_editing_level','bim_editing_level','sep_editing_level','squ_editing_level','lin_editing_level', 'bob_editing_level'], loc='lower left')
# plt.savefig(output_folder+'highly_conserved_recoding_sites_editing_levels_deviations.png')
# plt.close()
# 
# plt.gca().set_color_cycle(['yellow', 'blue','red','green','blue','black'])
# plt.plot(x,conserved_recoding_across_all['average_editing'],linewidth=0.5)
# plt.plot(x,conserved_recoding_across_all['lin_editing_level'],linewidth=0.5)
# plt.plot(x,conserved_recoding_across_all['bob_editing_level'],linewidth=0.5)
# plt.xlabel('site index')
# plt.ylabel('editing level')
# plt.legend(['average editing level','lin_editing_level','bob_editing_level'], loc='uper left')
# plt.savefig(output_folder+'editing_level_comparison.png')
# plt.close()
# 
# =============================================================================
