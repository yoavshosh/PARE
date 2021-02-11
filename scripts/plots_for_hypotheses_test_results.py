# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 20:12:28 2020

@author: shosh
"""

import os
import glob
import sys
import pandas as pd
import numpy as np
from functools import reduce
import itertools as it
from scipy import stats
from Bio import Phylo, SeqIO
from Bio.Seq import Seq

from math import floor, log10, sqrt, log
from decimal import Decimal
import gzip
import re

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import lines 
    import matplotlib.transforms as transforms
    import pylab
except ImportError:
    print('Could not Import data-viz libraries')

from hypotheses_test import *

try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3


#this dictionaty defines the labeles of each species in plots and results tables
animals_names_dict={'sep':'S.ofi',
                    'squ':'D.pea',
                    'bim':'O.bim',
                    'oct':'O.vul',
                    'bob':'Eup',
                    'lin':'S.lin',
                    'apl':'Apl',
                    'nau':'Nau',
                    'bim_oct':'O.vul-O.bim',
                    'sep_squ':'D.pea-S.ofi',
                    'lin_sep':'S.ofi-S.lin',
                    'lin_squ':'D.pea-S.lin',
                    'bob_lin_squ':'D.pea-S.lin-Eup',
                    'bob_lin':'Eup-S.lin',
                    'bob_lin_sep':'S.ofi-S.lin-Eup',
                    'bob_lin_sep_squ':'Decapodiformes',
                    'bim_bob_lin_oct_sep_squ':'Coleoids',
                    'N0':'Mollusca',
                    'N1':'C0',
                    'C':'C1',
                    'D':'D',
                    'B':'B',
                    'S':'S',
                    'O':'O',
                    'S0':'S0',
                    'S1':'S1'}


def read_fasta_to_df(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['id','sequence']
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        data.append((record.id.split('|')[1].rstrip().lstrip(),record.seq))
    
    df = pd.DataFrame(data = data, columns = col_names)
    return df


def find_codon(row, a, trinity_df):
    coding_loc=int(row[a+'_coding_location'])
    loc_in_codon=coding_loc%3
    sequence=trinity_df.loc[row[a+'_id'].split('|')[1],'sequence']
    codon=sequence[coding_loc-loc_in_codon:coding_loc+3-loc_in_codon]
    if list(codon).count('A')>1:
        codon=codon+';'+str(loc_in_codon)
    return str(codon)
        

def codons_dist_dict(sites_df,trinity_dict,conserved_animals,aa_instead=False):
    
    sites_df['take_site']=sites_df.apply(lambda row: all([int(row[a+'_edited']) for a in conserved_animals]), axis=1)
    relevant_sites = sites_df[sites_df['take_site']].copy()
    if aa_instead:
        relevant_sites['consensus_swap']=relevant_sites.apply(lambda row: max(set([row[a+'_original_aa']+row[a+'_AG_target_aa'] for a in conserved_animals]),
                                                                              key=[row[a+'_original_aa']+row[a+'_AG_target_aa'] for a in conserved_animals].count), axis=1)
    else:
        for a in conserved_animals:
            relevant_sites[a+'_codon']=relevant_sites.apply(lambda row: find_codon(row,a,trinity_dict[a]),axis=1)
        relevant_sites['consensus_swap']=relevant_sites.apply(lambda row: max(set([row[a+'_codon'] for a in conserved_animals]),
                                                                              key=[row[a+'_codon'] for a in conserved_animals].count), axis=1)
        
    swaps_dict={}
    for swap in sorted(set(relevant_sites['consensus_swap'].values)):
        swaps_dict.update({swap:len(relevant_sites[relevant_sites['consensus_swap']==swap])})    
    return swaps_dict
    

def mismatches_dist(outpath,codons_dist_dict,animals):
    
    def edit_codon(codon,loc):
        codon_lst=list(codon)
        codon_lst[int(loc)]='G'
        return ''.join(codon_lst)            
    
    codons_lst=[]
    for p in it.permutations((('A'),('C','T','G'),('A'),('C','T','G'),('A')),3):
        for c in it.product(*p):
            codon=''.join(c)
            for loc in [m.start() for m in re.finditer('A', codon)]:
                codons_lst.append(codon+';'+str(loc))
    codons_lst=set(codons_lst)
    
    data=[]
    if type(animals) is str:
        values_col = animals+'_events'
    else:
        values_col = '_'.join(animals)+'_events'
    columns = ('original_aa','target_aa','codon','nucl',values_col)
    for codon in codons_lst:
        codon_str=codon.split(';')[0]
        loc=int(codon.split(';')[1])
        try:
            events=codons_dist_dict[codon]
        except KeyError:
            events=0
        if len(codon_str)!=3:
            print(str(animals)+' '+codon+' is not a multiplete of 3')

        data.append((str(Seq(codon).translate()),
                     str(Seq(edit_codon(codon_str,loc)).translate()),
                     codon_str, loc, events))
        
    df = pd.DataFrame(data=data,columns=columns)
    sorted_df=df.sort_values(['original_aa','codon','nucl'],ascending=[True,True,True])
    return sorted_df

def plot_mm_dist(df,outpath,title):

        width=0.2
        sorted_df=df.sort_values(['original_aa','codon','nucl'],ascending=[True,True,True])
        sorted_df['editing_in_codon']=sorted_df.apply(lambda row: row['codon']+';'+str(row['nucl']+1)+'>'+row['target_aa'], axis=1)
        fig, ax = plt.subplots(figsize=(20, 10))
        swaps=list(sorted_df['editing_in_codon'].values)
        events = list(sorted_df[title+'_events'].values)
        values = [x/sum(events) for x in events]
        plt.bar(swaps,values,width=width)
        plt.xticks(np.arange(len(swaps)),swaps,rotation=90)   
        plt.ylabel('Mismatch fraction')
        aas=sorted(list(set(list(sorted_df['original_aa'].values))))
        ticks=ax.get_xticks()

        j=0
        k=0
        ticks_vals = list(enumerate(ticks))
        trans = ax.get_xaxis_transform()
        for i,t in ticks_vals:
            if i==j:
                aa=aas[k]
                ticks_n=len(sorted_df[sorted_df['original_aa']==aa])
                last_tick_for_aa=ticks_vals[i+ticks_n-1]
                plt.text((t+float(ticks_n-1)/2-0.15),-0.15,aa,size=10,transform=trans)    
                plt.plot([t-0.30,last_tick_for_aa[1]+0.25],[-0.13,-0.13],color = 'black',transform=trans,clip_on=False)        
                j+=ticks_n
                k+=1
            else:
                pass

        plt.savefig(outpath+title+'_codons_dist.jpg')
        plt.close()
        
        
def plot_mm_dist_for_multiple_animals(df,outpath,title,animals):

        sorted_df=df.sort_values(['original_aa','codon','nucl'],ascending=[True,True,True])
        sorted_df['editing_in_codon']=sorted_df.apply(lambda row: row['codon']+';'+str(row['nucl']+1)+'>'+row['target_aa'], axis=1)
        fig, ax = plt.subplots(figsize=(20, 10))
        
        
        colors = ['red','blue','yellow','grey','green','pink']
        
        bar_width=0.5
        swaps=list(sorted_df['editing_in_codon'].values)
        ticks_steps=4
        locs = np.arange(start=0,stop=len(swaps)*ticks_steps,step=ticks_steps)   
        for i,a in enumerate(animals):
            bars_locs=[x+i*bar_width for x in locs]
            events = list(sorted_df[a].values)
            values = [x/sum(events) for x in events]
            plt.bar(bars_locs,values,width=bar_width,color=colors[i],label=a)
        
        plt.xticks([x+bar_width*(len(animals)-1)/2 for x in locs],swaps,rotation=90)   
        plt.ylabel('Mismatch fraction',fontsize='20')
        plt.tick_params(axis='y', which='both', labelsize=17)
        aas=sorted(list(set(list(sorted_df['original_aa'].values))))
        ticks=ax.get_xticks()

        j=0
        k=0
        ticks_vals = list(enumerate(ticks))
        trans = ax.get_xaxis_transform()
        for i,t in ticks_vals:
            if i==j:
                aa=aas[k]
                ticks_n=len(sorted_df[sorted_df['original_aa']==aa])
                last_tick_for_aa=ticks_vals[i+ticks_n-1]
                plt.text((t+float(ticks_n-1)*ticks_steps/2-0.8),-0.15,aa,size=10,transform=trans)    
                plt.plot([t-1.2,last_tick_for_aa[1]+0.9],[-0.13,-0.13],color = 'black',transform=trans,clip_on=False)        
                j+=ticks_n
                k+=1
            else:
                pass
            
        plt.legend(loc=(0.85,0.6),prop={'size':20})
        plt.savefig(outpath+title+'.jpg')
        plt.close()
    


def collect_results_for_anacestral_state_hpm_and_plot_probs(path,animals,conserved_groups):
    
    def plot_rates_grouped_bars(path,df,name,title):
                
        df['Species']=df.apply(lambda row: animals_names_dict[row.name],axis=1)
        df.set_index('Species',inplace=True)
                
        df['syn_res_rate'] = df.apply(lambda row: float(row.syn_res_edited)/row.syn_res, axis=1)
        df['nonsyn_res_rate'] = df.apply(lambda row: float(row.nonsyn_res_edited)/row.nonsyn_res, axis=1)
        df['syn_div_rate'] = df.apply(lambda row: float(row.syn_div_edited)/row.syn_div, axis=1)
        df['nonsyn_div_rate'] = df.apply(lambda row: float(row.nonsyn_div_edited)/row.nonsyn_div, axis=1)
        df['syn_res_err'] = df.apply(lambda row: np.sqrt(row['syn_res_rate']*(1-row['syn_res_rate'])/row.syn_res), axis=1)
        df['nonsyn_res_err'] = df.apply(lambda row: np.sqrt(row['nonsyn_res_rate']*(1-row['nonsyn_res_rate'])/row.nonsyn_res), axis=1)
        df['syn_div_err'] = df.apply(lambda row: np.sqrt(row['syn_div_rate']*(1-row['syn_div_rate'])/row.syn_div), axis=1)
        df['nonsyn_div_err'] = df.apply(lambda row: np.sqrt(row['nonsyn_div_rate']*(1-row['nonsyn_div_rate'])/row.nonsyn_div), axis=1)
        df['res_div_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_res_edited,row.nonsyn_res],[row.nonsyn_div_edited,row['nonsyn_div']]])[1], axis=1)        
        df['res_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_res_edited,row.nonsyn_res],[row.syn_res_edited,row['syn_res']]])[1], axis=1)        
        df['div_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_div_edited,row['nonsyn_div']],[row.syn_div_edited,row['syn_div']]])[1], axis=1)         
        max_y = max(list(df['synonymous'])+list(df['restorative'])+list(df['diversifying']))        
        pdax = pd.concat([df.synonymous.rename('synonymous'), df.restorative.rename('restorative'), df.diversifying.rename('diversifying')], axis=1).plot(kind='bar', yerr=(df['syn_err'].values,df['res_err'].values,df['div_err'].values), error_kw=dict(lw=1, capsize=2, capthick=1), color=['grey','deepskyblue','red'])
        ticks = pdax.get_xticks()
        labels = pdax.get_xticklabels()
        plt.xticks(rotation=0)
        plt.title(title+'\n\n', fontsize=17)
        plt.ylabel('Fraction of sites edited', fontsize=17)
        plt.xlabel('Species', fontsize=17)
        plt.tick_params(axis='x', which='both', length=0, labelsize=13)
        plt.tick_params(axis='y', which='both', labelsize=15)
        plt.ylim([0,max_y*1.2])
        plt.gcf().set_size_inches(8,7)
        pdax.legend(ncol=3,fontsize=13,bbox_to_anchor=(0.95, 1.1))
        for i,t in enumerate(ticks):
            plt.plot([t,t+0.2],[max_y*1.08,max_y*1.08],color = 'black')
            plt.plot([t,t],[max_y*1.08,max_y*1.07],color = 'black')
            plt.plot([t+0.2,t+0.2],[max_y*1.08,max_y*1.07],color = 'black')
            ps=pval_str(df.loc[labels[i].get_text(),'all_sites_pval'])
            if ps=='ns':
                align=0.03
            else:
                align=-0.15
            plt.text(t+align,max_y*1.1,pval_str(df.loc[labels[i].get_text(),'all_sites_pval']),size=10)
        plt.savefig(path+name+'_types_rates.png')
        plt.close()
    
        df['syn_res_species_specific_rate'] = df.apply(lambda row: float(row.syn_res_species_specific_edited)/row.syn_res, axis=1)
        df['nonsyn_res_species_specific_rate'] = df.apply(lambda row: float(row.nonsyn_res_species_specific_edited)/row.nonsyn_res, axis=1)
        df['syn_div_species_specific_rate'] = df.apply(lambda row: float(row.syn_div_species_specific_edited)/row.syn_div, axis=1)
        df['nonsyn_div_species_specific_rate'] = df.apply(lambda row: float(row.nonsyn_div_species_specific_edited)/row.nonsyn_div, axis=1)
        df['syn_res_species_specific_err'] = df.apply(lambda row: np.sqrt(row['syn_res_species_specific_rate']*(1-row['syn_res_species_specific_rate'])/row.syn_res), axis=1)
        df['nonsyn_res_species_specific_err'] = df.apply(lambda row: np.sqrt(row['nonsyn_res_species_specific_rate']*(1-row['nonsyn_res_species_specific_rate'])/row.nonsyn_res), axis=1)
        df['syn_div_species_specific_err'] = df.apply(lambda row: np.sqrt(row['syn_div_species_specific_rate']*(1-row['syn_div_species_specific_rate'])/row.syn_div), axis=1)
        df['nonsyn_div_species_specific_err'] = df.apply(lambda row: np.sqrt(row['nonsyn_div_species_specific_rate']*(1-row['nonsyn_div_species_specific_rate'])/row.nonsyn_div), axis=1)
        df['res_div_species_specific_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_res_edited,row.nonsyn_res],[row.nonsyn_div_edited,row['nonsyn_div']]])[1], axis=1)        
        df['res_species_specific_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_res_edited,row.nonsyn_res],[row.syn_res_edited,row['syn_res']]])[1], axis=1)        
        df['div_species_specific_all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.nonsyn_div_edited,row['nonsyn_div']],[row.syn_div_edited,row['syn_div']]])[1], axis=1)        
        max_y = max(list(df['synonymous'])+list(df['restorative'])+list(df['diversifying']))
        pdax = pd.concat([df.synonymous_species_specific.rename('synonymous'), df.restorative_species_specific.rename('restorative'), df.diversifying_species_specific.rename('diversifying')], axis=1).plot(kind='bar',yerr=(df['syn_ss_err'].values,df['res_ss_err'].values,df['div_ss_err'].values), error_kw=dict(lw=1, capsize=2, capthick=1),color=['grey','deepskyblue','red'])           
        ticks = pdax.get_xticks()
        labels = pdax.get_xticklabels()
        plt.xticks(rotation=0)
        plt.title(title+'\n\n', fontsize=17)
        plt.ylabel('Fraction of sites edited', fontsize=17)
        plt.xlabel('Species', fontsize=17)
        plt.tick_params(axis='x', which='both', length=0, labelsize=13)
        plt.tick_params(axis='y', which='both', labelsize=15)   
        plt.ylim([0,max_y*1.2])
        plt.gcf().set_size_inches(8,7)
        pdax.legend(ncol=3,fontsize=13,bbox_to_anchor=(0.95, 1.1))
        for i,t in enumerate(ticks):
            plt.plot([t,t+0.2],[max_y*1.08,max_y*1.08],color = 'black')
            plt.plot([t,t],[max_y*1.08,max_y*1.07],color = 'black')
            plt.plot([t+0.2,t+0.2],[max_y*1.08,max_y*1.07],color = 'black')
            ps=pval_str(df.loc[labels[i].get_text(),'all_sites_pval'])
            if ps=='ns':
                align=0
            else:
                align=-0.15
            plt.text(t+align,max_y*1.1,pval_str(df.loc[labels[i].get_text(),'all_sites_pval']),size=10)
        plt.savefig(path+'species_specific_'+name+'_types_rates.png')
        plt.close()
        return df
        
    def plot_editing_levels_distributions_by_types(path,el_by_types_dict,animals,name,ylim):
        fig, axes = plt.subplots(nrows=1, ncols=len(animals), figsize=(len(animals)*3, 7))
        for i,a in enumerate(animals):    
            colors = ['grey','deepskyblue','red']
            data = [el_by_types_dict[a+'_syn'].values,el_by_types_dict[a+'_res'].values,el_by_types_dict[a+'_div'].values]
            pval = stats.mannwhitneyu(el_by_types_dict[a+'_res'].values,el_by_types_dict[a+'_div'].values)[1]
            pval_s = pval_str(pval)
            bplot = axes[i].boxplot(data,showfliers=False,patch_artist=True,labels=['syn','res','div'])
            axes[i].set_title(animals_names_dict[a], fontsize=17,loc='Center')
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)
            axes[i].set_ylim(0,ylim)
            medians = [np.median(x) for x in data]
            pos = np.arange(3) + 1
            upper_labels = [str(np.round(s, 4)) for s in medians]
            for tick, label in zip(range(3), axes[i].get_xticklabels()):
                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round_siginficant(float(u))) for j, u in enumerate(upper_labels)])
                
                label.set_rotation(90)
                axes[i].tick_params(axis='x', which='both', length=0, labelsize=15)
            if i:
                axes[i].set_yticks([])
            else:
                axes[i].set_yticks([0,0.1,0.2,0.3,0.4])
                axes[i].set_yticklabels([0,0.1,0.2,0.3,0.4],rotation=90)
                axes[i].tick_params(axis='y', which='both', length=0, labelsize=15)
            if not i:
                axes[i].set_ylabel('Editing level', fontsize=17)
            axes[i].axhline(y=ylim-0.05,xmin=0.45,xmax=0.9,c="black",linewidth=0.5)
            
            align=-0.4
            axes[i].text(2.35+align,ylim-0.045,pval_s,size=15)

        fig.savefig(path+name+'_sites_editing_level_dist_by_types.png')
        plt.close()
        
        fig, axes = plt.subplots(nrows=1, ncols=len(animals), figsize=(len(animals)*3, 7))
        for i,a in enumerate(animals):    
            colors = ['grey','deepskyblue','red']
            data = [el_by_types_dict[a+'_specific_syn'].values,el_by_types_dict[a+'_specific_res'].values,el_by_types_dict[a+'_specific_div'].values]
            pval = stats.mannwhitneyu(el_by_types_dict[a+'_specific_res'].values,el_by_types_dict[a+'_specific_div'].values)[1]
            pval_s = pval_str(pval)
            bplot = axes[i].boxplot(data,showfliers=False,patch_artist=True,labels=['syn','res','div'])
            axes[i].set_title(animals_names_dict[a], fontsize=17,loc='Center')
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)
            axes[i].set_ylim(0,ylim)
            medians = [np.median(x) for x in data]
            pos = np.arange(3) + 1
            upper_labels = [str(np.round(s, 4)) for s in medians]
            for tick, label in zip(range(3), axes[i].get_xticklabels()):
                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round_siginficant(float(u))) for j, u in enumerate(upper_labels)])
                
                label.set_rotation(90)
                axes[i].tick_params(axis='x', which='both', length=0, labelsize=15)
            if i:
                axes[i].set_yticks([])
            else:
                axes[i].set_yticks([0,0.1,0.2,0.3,0.4])
                axes[i].set_yticklabels([0,0.1,0.2,0.3,0.4],rotation=90)
                axes[i].tick_params(axis='y', which='both', length=0, labelsize=15)
            if not i:
                axes[i].set_ylabel('Editing level', fontsize=17)
            axes[i].axhline(y=ylim-0.05,xmin=0.45,xmax=0.9,c="black",linewidth=0.5)
            align=-0.4
            axes[i].text(2.35+align,ylim-0.045,pval_s,size=15)
            
        fig.savefig(path+name+'_species_specific_sites_editing_level_dist_by_types.png')
        plt.close()
        
    #editing types count - just convert to xlsx
    cnt = pd.read_csv(path+'editing_types_count',sep='\t',index_col=0)
    cnt.to_excel(path+'editing_types_count.xlsx')

    #editing types rates
    # rates = pd.read_csv(path+'editing_types_rates',sep='\t',index_col=0)
    # rates.to_excel(path+'editing_types_rates.xlsx')
    # plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites')
    # plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites')
    
        # editing levels by types
    # el_by_types_dict={}
    # for f in glob.glob(path+'editing_levels_by_types_*'):
        # name=f.split('\\')[-1].replace('editing_levels_by_types_','')
        # s = pd.read_csv(f, sep='\t', index_col=0, squeeze=True)
        # el_by_types_dict.update({name:s})
    # plot_editing_levels_distributions_by_types(path,el_by_types_dict,animals,'all',0.4)
    # plot_editing_levels_distributions_by_types(path,el_by_types_dict,conserved_groups,'conserved',1)
    
    #editing types rates
    rates = pd.read_csv(path+'ancestral_editing_types_rates_strong',sep='\t',index_col=0)
    rates_claced = plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites_strong','Strong sites')
    rates_claced.to_excel(path+'ancestral_editing_types_rates_strong.xlsx')
    rates_claced = plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites_strong','Strong sites')
    rates_claced.to_excel(path+'ancestral_conserved_editing_types_rates_strong.xlsx')
    # df.to_excel(path+'ancestral_conserved_editing_types_rates_strong.xlsx')
    
    #editing types rates
    rates = pd.read_csv(path+'editing_types_rates_weak',sep='\t',index_col=0)
    rates_claced = plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites_weak','Weak sites')
    rates_claced.to_excel(path+'editing_types_rates_weak.xlsx')
    rates_claced = plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites_weak','Weak sites')
    rates_claced.to_excel(path+'conserved_editing_types_rates_weak.xlsx')


def read_editing_sites_tbl(editing_sites_file, mm_type = 'AG'):
    """
    read the editing sites tabel
    """
    col_names = ['id', 'protein', 'location', 'mm_type', 'DNA_A', 'DNA_T', 'DNA_G', 'DNA_C',
                 'RNA_A', 'RNA_T', 'RNA_G', 'RNA_C', 'Trinity', 'RNA_coverage', 'DNA_coverage', 
                 'p_val', 'AA_before', 'AA_after', 'type', 'protein_length', 'editing_level', 'strand']
    sites_df = pd.read_csv(editing_sites_file, sep = '\t', names = col_names, index_col = None)
    sites_df = sites_df[sites_df['mm_type'] == mm_type]
    
    return sites_df


def read_trinity_mrna_files(trinity_file):
    """
    read relevat information from transcriptome file
    """
    data = []
    col_names = ['id','protein','protein_name','orfs_start','orfs_end','strand','sequence']
    
    for record in SeqIO.parse(open(trinity_file, "r"), "fasta"):
        rec_data = record.description.split('\t')
        if rec_data[0][-1] == ' ':  #some fastafiles have spaces after each id, so fixing it here.
            rec_data[0] = rec_data[0].replace(' ','')
        protein = rec_data[-1].split('|')[2].split(' ')[0]   #reading proteing from description assuming it was added to header using the transcriptome built pipeline we have for trinity
        rec_data = (rec_data[0],protein,protein.split('_')[0],int(rec_data[2]),int(rec_data[4]),rec_data[6],record.seq)
        data.append(rec_data)
    
    df = pd.DataFrame(data = data, columns = col_names)
    
    return df


def draw_tree(tree_str,out_path):
    handle = StringIO(tree_str)
    tree = Phylo.read(handle, "newick")
    tree.ladderize()
    fig = plt.figure(figsize=(20, 10), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    #Phylo.draw(tree, axes=axes,do_show=False,branch_labels=lambda c: c.branch_length)
    Phylo.draw(tree, axes=axes,do_show=True)
    pylab.savefig(out_path+'tree.png',format='png', bbox_inches='tight')
    

def pval_str(p,s=2):

    if p==0:
        return '0 *'
    elif pd.isnull(p):
        return ''
    elif p<0.05:
        # print(p)
        p=round(float(p),s-int(floor(log10(abs(float(p)))))-1)
        return '%.1e' % Decimal(p) + ' *'
    else:
        return '%.1e' % Decimal(p)


def round_siginficant(n,s=2):
    r=round(float(n),s-int(floor(log10(abs(float(n)))))-1)
    r_str=str(r).replace('.','')
    while r_str[0]=='0':
        r_str=r_str[1:]
    zeroes_to_add = s-len(r_str)
    
    new_r_str=str(r)
    for i in range(zeroes_to_add):
        new_r_str=new_r_str+'0'    
    
    return new_r_str


def pval_stars(p):
    if p>0.05:
        pvals_str='ns'
    elif p<0.0001:
        pvals_str='****'
    elif p<0.001:
        pvals_str='***'
    elif p<0.01:
        pvals_str='**'
    elif p<0.05:
        pvals_str='*'
    return pvals_str

    
def collect_results_for_hpm_and_plot_probs(path,animals,conserved_groups=None):
        
    def plot_rates_grouped_bars(path,df,name,title):
                
        df['Species']=df.apply(lambda row: animals_names_dict[row.name],axis=1)
        df.set_index('Species',inplace=True)
        df['synonymous'] = df.apply(lambda row: float(row.syn_edited)/row.syn, axis=1)
        df['restorative'] = df.apply(lambda row: float(row.res_edited)/row.res, axis=1)
        df['diversifying'] = df.apply(lambda row: float(row.div_edited)/row['div'], axis=1)
        df['syn_err'] = df.apply(lambda row: np.sqrt(row['synonymous']*(1-row['synonymous'])/row.syn), axis=1)
        df['res_err'] = df.apply(lambda row: np.sqrt(row['restorative']*(1-row['restorative'])/row.res), axis=1)
        df['div_err'] = df.apply(lambda row: np.sqrt(row['diversifying']*(1-row['diversifying'])/row['div']), axis=1)
        df['all_sites_pval'] = df.apply(lambda row: stats.fisher_exact([[row.res_edited,row.res],[row.div_edited,row['div']]])[1], axis=1)        
        max_y = max(list(df['synonymous'])+list(df['restorative'])+list(df['diversifying']))        
        pdax = pd.concat([df.synonymous.rename('synonymous'), df.restorative.rename('restorative'), df.diversifying.rename('diversifying')], axis=1).plot(kind='bar', yerr=(df['syn_err'].values,df['res_err'].values,df['div_err'].values), error_kw=dict(lw=1, capsize=2, capthick=1), color=['grey','deepskyblue','red'])
        ticks = pdax.get_xticks()
        labels = pdax.get_xticklabels()
        plt.xticks(rotation=0)
        plt.title(title+'\n\n', fontsize=17)
        plt.ylabel('Fraction of sites edited', fontsize=17)
        plt.xlabel('Species', fontsize=17)
        plt.tick_params(axis='x', which='both', length=0, labelsize=13)
        plt.tick_params(axis='y', which='both', labelsize=15)
        plt.ylim([0,max_y*1.2])
        plt.gcf().set_size_inches(8,7)
        pdax.legend(ncol=3,fontsize=13,bbox_to_anchor=(0.95, 1.1))
        for i,t in enumerate(ticks):
            plt.plot([t,t+0.2],[max_y*1.08,max_y*1.08],color = 'black')
            plt.plot([t,t],[max_y*1.08,max_y*1.07],color = 'black')
            plt.plot([t+0.2,t+0.2],[max_y*1.08,max_y*1.07],color = 'black')
            ps=pval_str(df.loc[labels[i].get_text(),'all_sites_pval'])
            if ps=='ns':
                align=0.03
            else:
                align=-0.15
            plt.text(t+align,max_y*1.1,pval_str(df.loc[labels[i].get_text(),'all_sites_pval']),size=10)
        plt.savefig(path+name+'_types_rates.png')
        plt.close()
        
        df['synonymous_species_specific'] = df.apply(lambda row: float(row.specie_specific_syn_edited)/row.syn, axis=1)
        df['restorative_species_specific'] = df.apply(lambda row: float(row.specie_specific_res_edited)/row.res, axis=1)
        df['diversifying_species_specific'] = df.apply(lambda row: float(row.specie_specific_div_edited)/row['div'], axis=1)
        df['syn_ss_err'] = df.apply(lambda row: np.sqrt(row['synonymous']*(1-row['synonymous'])/row.syn), axis=1)
        df['res_ss_err'] = df.apply(lambda row: np.sqrt(row['restorative']*(1-row['restorative'])/row.res), axis=1)
        df['div_ss_err'] = df.apply(lambda row: np.sqrt(row['diversifying']*(1-row['diversifying'])/row['div']), axis=1)
        df['species_specific_pval'] = df.apply(lambda row: stats.fisher_exact([[row.res_edited,row.res],[row.div_edited,row['div']]])[1], axis=1)
        max_y = max(list(df['synonymous'])+list(df['restorative'])+list(df['diversifying']))
        pdax = pd.concat([df.synonymous_species_specific.rename('synonymous'), df.restorative_species_specific.rename('restorative'), df.diversifying_species_specific.rename('diversifying')], axis=1).plot(kind='bar',yerr=(df['syn_ss_err'].values,df['res_ss_err'].values,df['div_ss_err'].values), error_kw=dict(lw=1, capsize=2, capthick=1),color=['grey','deepskyblue','red'])           
        ticks = pdax.get_xticks()
        labels = pdax.get_xticklabels()
        plt.xticks(rotation=0)
        plt.title(title+'\n\n', fontsize=17)
        plt.ylabel('Fraction of sites edited', fontsize=17)
        plt.xlabel('Species', fontsize=17)
        plt.tick_params(axis='x', which='both', length=0, labelsize=13)
        plt.tick_params(axis='y', which='both', labelsize=15)   
        plt.ylim([0,max_y*1.2])
        plt.gcf().set_size_inches(8,7)
        pdax.legend(ncol=3,fontsize=13,bbox_to_anchor=(0.95, 1.1))
        for i,t in enumerate(ticks):
            plt.plot([t,t+0.2],[max_y*1.08,max_y*1.08],color = 'black')
            plt.plot([t,t],[max_y*1.08,max_y*1.07],color = 'black')
            plt.plot([t+0.2,t+0.2],[max_y*1.08,max_y*1.07],color = 'black')
            ps=pval_str(df.loc[labels[i].get_text(),'all_sites_pval'])
            if ps=='ns':
                align=0
            else:
                align=-0.15
            plt.text(t+align,max_y*1.1,pval_str(df.loc[labels[i].get_text(),'all_sites_pval']),size=10)
        plt.savefig(path+'species_specific_'+name+'_types_rates.png')
        plt.close()
        return df
        
    def plot_editing_levels_distributions_by_types(path,el_by_types_dict,animals,name,ylim):
        fig, axes = plt.subplots(nrows=1, ncols=len(animals), figsize=(len(animals)*3, 7))
        for i,a in enumerate(animals):    
            colors = ['grey','deepskyblue','red']
            data = [el_by_types_dict[a+'_syn'].values,el_by_types_dict[a+'_res'].values,el_by_types_dict[a+'_div'].values]
            pval = stats.mannwhitneyu(el_by_types_dict[a+'_res'].values,el_by_types_dict[a+'_div'].values)[1]
            pval_s = pval_str(pval)
            bplot = axes[i].boxplot(data,showfliers=False,patch_artist=True,labels=['syn','res','div'],medianprops=dict(color='black'))
            axes[i].set_title(animals_names_dict[a], fontsize=17,loc='Center')
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)
            axes[i].set_ylim(0,ylim)
            # print(a)
            medians = [np.median(x) for x in data]
            pos = np.arange(3) + 1
            upper_labels = [str(np.round(s, 4)) for s in medians]
            for tick, label in zip(range(3), axes[i].get_xticklabels()):
                # axes[i].text(pos[tick], ylim-0.05, upper_labels[tick] ,horizontalalignment='center', verticalalignment='baseline', color='black') 
                # axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round(float(u),3-int(floor(log10(abs(float(u)))))-1)) for j, u in enumerate(upper_labels)])
                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round_siginficant(float(u))) for j, u in enumerate(upper_labels)])
                # axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(u) for j, u in enumerate(upper_labels)])
                
                label.set_rotation(90)
                axes[i].tick_params(axis='x', which='both', length=0, labelsize=15)
            if i:
                axes[i].set_yticks([])
            else:
                axes[i].set_yticks([0,0.1,0.2,0.3,0.4])
                axes[i].set_yticklabels([0,0.1,0.2,0.3,0.4],rotation=90)
                axes[i].tick_params(axis='y', which='both', length=0, labelsize=15)
            if not i:
                axes[i].set_ylabel('Editing level', fontsize=17)
            axes[i].axhline(y=ylim-0.05,xmin=0.45,xmax=0.9,c="black",linewidth=0.5)
            
            align=-0.4
            
            axes[i].text(2.35+align,ylim-0.045,pval_s,size=15)

        fig.savefig(path+name+'_sites_editing_level_dist_by_types.png')
        plt.close()

        fig, axes = plt.subplots(nrows=1, ncols=len(animals), figsize=(len(animals)*3, 7))
        for i,a in enumerate(animals):    
            colors = ['grey','deepskyblue','red']
            data = [el_by_types_dict[a+'_specific_syn'].values,el_by_types_dict[a+'_specific_res'].values,el_by_types_dict[a+'_specific_div'].values]
            pval = stats.mannwhitneyu(el_by_types_dict[a+'_specific_res'].values,el_by_types_dict[a+'_specific_div'].values)[1]
            pval_s = pval_str(pval)
            bplot = axes[i].boxplot(data,showfliers=False,patch_artist=True,labels=['syn','res','div'])
            axes[i].set_title(animals_names_dict[a], fontsize=17,loc='Center')
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)
            axes[i].set_ylim(0,ylim)
            medians = [np.median(x) for x in data]
            pos = np.arange(3) + 1
            upper_labels = [str(np.round(s, 4)) for s in medians]
            for tick, label in zip(range(3), axes[i].get_xticklabels()):
#                axes[i].text(pos[tick], ylim-0.05, upper_labels[tick] ,horizontalalignment='center', verticalalignment='baseline', color='black') 
#                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round(float(u),3-int(floor(log10(abs(float(u)))))-1)) for j, u in enumerate(upper_labels)])
                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(round_siginficant(float(u))) for j, u in enumerate(upper_labels)])
#                axes[i].set_xticklabels([['syn\n','res\n','div\n'][j]+str(u) for j, u in enumerate(upper_labels)])
                
                label.set_rotation(90)
                axes[i].tick_params(axis='x', which='both', length=0, labelsize=15)
            if i:
                axes[i].set_yticks([])
            else:
                axes[i].set_yticks([0,0.1,0.2,0.3,0.4])
                axes[i].set_yticklabels([0,0.1,0.2,0.3,0.4],rotation=90)
                axes[i].tick_params(axis='y', which='both', length=0, labelsize=15)
            if not i:
                axes[i].set_ylabel('Editing level', fontsize=17)
            axes[i].axhline(y=ylim-0.05,xmin=0.45,xmax=0.9,c="black",linewidth=0.5)
            
            align=-0.4

            
            axes[i].text(2.35+align,ylim-0.045,pval_s,size=15)
#            plt.xticks(rotation=90)

        fig.savefig(path+name+'_species_specific_sites_editing_level_dist_by_types.png')
        plt.close()
        

#    #editing types rates
#    rates = pd.read_csv(path+'editing_types_rates',sep='\t',index_col=0)
#    rates.to_excel(path+'editing_types_rates.xlsx')
#    plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites')
#    plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites')
    
    rates_dict = {}
    #editing types rates
    rates = pd.read_csv(path+'editing_types_rates_strong',sep='\t',index_col=0)
    rates_claced = plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites_strong','Strong sites')
    rates_claced.to_excel(path+'editing_types_rates_strong.xlsx')
    rates_dict.update({'strong':rates_claced})
    if conserved_groups is not None:
        rates_claced = plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites_strong','Strong sites')
        rates_claced.to_excel(path+'conserved_editing_types_rates_strong.xlsx')
    
    
        #editing types rates
    rates = pd.read_csv(path+'editing_types_rates_weak',sep='\t',index_col=0)
    rates_claced = plot_rates_grouped_bars(path,rates.loc[animals,:],'all_sites_weak','Weak sites')
    rates_claced.to_excel(path+'editing_types_rates_weak.xlsx')
    rates_dict.update({'weak':rates_claced})
    if conserved_groups is not None:
        rates_claced = plot_rates_grouped_bars(path,rates.loc[conserved_groups,:],'conserved_sites_weak','Strong sites')
        rates_claced.to_excel(path+'conserved_editing_types_rates_weak.xlsx')
    

    #editing levels by types
    el_by_types_dict={}
    for f in glob.glob(path+'editing_levels_by_type*'):
        name=f.split('\\')[-1].replace('editing_levels_by_types_','')
        try:
            s = pd.read_csv(f, sep='\t', index_col=0, header=None, squeeze=True)
        except pd.errors.EmptyDataError:
            s = None
        el_by_types_dict.update({name:s})
    plot_editing_levels_distributions_by_types(path,el_by_types_dict,animals,'all',0.4)
    if conserved_groups is not None:
        plot_editing_levels_distributions_by_types(path,el_by_types_dict,conserved_groups,'conserved',1)    
    
    
    return rates_dict, el_by_types_dict


def collect_results_for_general_model_and_plot_probs(path, recalc=True, iden_list=[0.3],iden_range_list=[10],editing_level_method_list=['average'],strong_levels=[[0.05,1],[0.1,1],[0.15,1],[0.2,1]], weak_levels=[[0,0.01],[0,0.02],[0,0.05],[0,0.1]],models=['adaptive']):
    
    def convert_dfs_csv_to_xlsx(*partial_paths_list, file_name='analysis'):
    # def convert_dfs_csv_to_xlsx(*partial_paths_list, file_name='analysis'):
        models_dict = {}
        files = []
        for p in partial_paths_list:
            files+=glob.glob('/'.join(p.split('/')[:-1])+'/'+p.split('/')[-1]+'*')
        for f in files:
            df = pd.read_csv(f, sep='\t', error_bad_lines=False, index_col=0, dtype='unicode')
            df = df.apply(pd.to_numeric, errors='coerce').fillna(df)
            models_dict.update({f.split('/')[-1].split('\\')[-1]:df})
        return models_dict

    def convert_models_dfs_back_to_series(models_dict,new_models_dict=None):        
        if new_models_dict is None:
            new_models_dict = {}
        else:
            assert(type(new_models_dict)==dict),"new_models_dict have to be a dictionary"
            
        for k, df in models_dict.items():      
            
            for col in df.columns:
                s = df[col]       
                index = [''.join([j for j in i if not j.isdigit()]).replace('-','').replace('.','').replace('liberal_','').replace('strict_','') for i in s.index]                
                name = k+'_'+s['combined_editing_level_method']+'_strong_lower'+str(s['strong_levels_lower_limit'])+'_weak_upper'+str(s['weak_levels_upper_limit'])+'_adapR'+str(s['adaptive_rate'])
                       
                new_s = pd.Series(data = s.values, index=index, name=name)
                new_s = new_s.reindex(index = ['ancestor','intermediate','leaf','optimized_adaptive_rate','optimized_percentile_of_observed_mutation',
                                               'combined_editing_level_method','strong_levels_lower_limit','strong_levels_upper_limit','weak_levels_lower_limit','weak_levels_upper_limit',
                                               'adaptive_rate','intermediate_syn_a','intermediate_non_syn_a','syn_ag_mut_in_leaf','non_syn_ag_mut_in_leaf',
                                               'syn_ag_mut_rate','non_syn_ag_mut_rate',
                                               'groups_of_edited_leaves','non_edited_leaves','weak_syn_sites','weak_nonsyn_sites',
                                               'nonsyn_sites_depletion_factor','strong_syn_sites','strong_nonsyn_sites','adaptive_sites',
                                               'expected_hp_strong_nonsyn_sites','strong_nonsyn_sites_excess',
                                               'expected_nonsyn_EG_mutations_from_adaptive','expected_nonsyn_EG_mutations_from_expected_hp','expected_nonsyn_EG_mutations_from_excess',
                                               'total_expected_EG_mutations','std','observed_strong_nonsyn_eg_mutations','p','strong_nonsyn_eg_mutations_rate',
                                               'observed_strong_syn_eg_mutations','strong_syn_eg_mutations_rate','excess_of_unmutated_strong_nonsyn_sites'])
                new_models_dict.update({name:new_s})
        return new_models_dict
    
    
    def plot_p_values_matrix(path,name,results_df,intermediates_to_group = ['S','B']):
        
        title = name 
        data = []
        
        for i, row in results_df.iterrows():
            if row['intermediate'] in intermediates_to_group:
                model_name = row['combined_editing_level_method']+'_'+row['ancestor']+'_to_'+'/'.join(intermediates_to_group)+'_to_leaf'
            else:
                model_name = row['combined_editing_level_method']+'_'+row['ancestor']+'_to_'+row['intermediate']+'_to_leaf'
            
            levels_and_leaf = str(row['weak_levels_upper_limit'])+'_'+str(row['strong_levels_lower_limit'])+'_'+row['leaf']+'+adapR'+str(row['adaptive_rate'])
            p=row['p']
            if p == 'lam < 0':
                p = np.nan
            data+=[(model_name,levels_and_leaf,p,)]
        
        df = pd.DataFrame(data=data, columns=['branches','levels_bounds_and_leaves','p'])
        df['p']=df.apply(pd.to_numeric, errors='coerce')['p']
        pivot_df = df.pivot("branches", "levels_bounds_and_leaves", "p")
        plt.figure(figsize=(45,12))
        ax = sns.heatmap(pivot_df,annot=True,vmin=0.0, vmax=0.1,cmap='coolwarm',linewidths=.5)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        plt.savefig(path+title+'.png')
        plt.close()
    
    def restor_and_plot_distributions(path, df, name, confidence_level=1, n=1000000, recalc=True):
        
        def calc_expected_mutations_distribution(expected_mut,excess_mut,adaptive_mut,
                                                 syn_a,non_syn_a,syn_ag_mut,non_syn_ag_mut,weak_syn_sites,weak_non_syn_sites,
                                                 strong_syn_sites,strong_hp_nonsyn_sites,adaptive_sites,confidence_level,
                                                 n,fix_depletion,calc_expected_hp_sites_mutations,calc_adaptive_sites_mutations):
            
            syn_ag_mut_rate = float(syn_ag_mut)/syn_a
            non_syn_ag_mut_rate = float(non_syn_ag_mut)/non_syn_a
            z = 1-(1-confidence_level)/2
            p_weak_nss=float(weak_non_syn_sites)/non_syn_a
            p_weak_ss=float(weak_syn_sites)/syn_a
            p_strong_ss=float(strong_syn_sites)/syn_a
            if fix_depletion is not None:
                depletion_factor = fix_depletion
            
            expected_non_syn_sites_dist = []
            for i in range(n):
                if fix_depletion is None:
                    p_weak_nss_rand = np.random.normal(p_weak_nss,z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a))
                    while p_weak_nss_rand<0 or p_weak_nss_rand>1:
                        p_weak_nss_rand = np.random.normal(p_weak_nss,z*np.sqrt(p_weak_nss*(1-p_weak_nss)/non_syn_a))
                    p_weak_ss_rand = np.random.normal(p_weak_ss,z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a))
                    while p_weak_ss_rand<0 or p_weak_nss_rand>1:
                        p_weak_ss_rand = np.random.normal(p_weak_ss,z*np.sqrt(p_weak_ss*(1-p_weak_ss)/syn_a)) 
                    depletion_factor=(p_weak_nss_rand/p_weak_ss_rand)
                p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
                while p_strong_ss_rand<0 or p_strong_ss_rand>1:
                    p_strong_ss_rand = np.random.normal(p_strong_ss,z*np.sqrt(p_strong_ss*(1-p_strong_ss)/syn_a))
                
                lam = depletion_factor*p_strong_ss_rand*non_syn_a
                if lam<strong_hp_nonsyn_sites:
                    expected_non_syn_sites_dist.append(np.random.poisson(lam))
                else:
                    expected_non_syn_sites_dist.append(strong_hp_nonsyn_sites)
                
            for i in range(n):
                
                #calculate mean of expected mut from expected nonsyn sites, if calculation flag is on, and a mean is within feasible range. if not assume 0 mutations
                p_non_syn_ag_mut_rate_rand = np.random.normal(non_syn_ag_mut_rate,z*np.sqrt(non_syn_ag_mut_rate*(1-non_syn_ag_mut_rate)/non_syn_a))
                expected_eg_mut_from_expected_nss = expected_non_syn_sites_dist[i]*p_non_syn_ag_mut_rate_rand
                if expected_eg_mut_from_expected_nss>0.0 and calc_expected_hp_sites_mutations:
                    expected_eg_mut_from_expected_nss_rand = np.random.poisson(expected_eg_mut_from_expected_nss)
                else:
                    expected_eg_mut_from_expected_nss_rand=0.0
                
                #calculate mean of expected mut from excess of nonsyn sites, if mean is within feasible range. if not assume 0 mutations
                expected_eg_mut_from_excess = (strong_hp_nonsyn_sites-expected_non_syn_sites_dist[i])*syn_ag_mut_rate
                if expected_eg_mut_from_excess>0:
                    expected_eg_mut_from_excess_rand=np.random.poisson(expected_eg_mut_from_excess)
                else:
                    expected_eg_mut_from_excess_rand=0.0
                
                #calculate mean of expected mut from adaptive sites, if calculation flag is on, and a mean is within feasible range. if not assume 0 mutations
                expected_adptive_sites_mutations = adaptive_sites*p_non_syn_ag_mut_rate_rand
                if expected_adptive_sites_mutations>0 and calc_adaptive_sites_mutations:
                    expected_eg_mut_from_adaptive_rand=np.random.poisson(expected_adptive_sites_mutations)
                else:
                    expected_eg_mut_from_adaptive_rand=0.0
                    
                expected_mut.append(expected_eg_mut_from_expected_nss_rand)
                excess_mut.append(expected_eg_mut_from_excess_rand)
                adaptive_mut.append(expected_eg_mut_from_adaptive_rand)
                
    

        def update_table_and_plot_values(row,name,outpath,n,confidence_level,recalc):
            
            if eval(row['optimized_adaptive_rate']) or any([row['strong_levels_lower_limit']!=0.1, row['weak_levels_upper_limit']!=0.05]):
                recalc=False
            
            if recalc:
                fig_name=name+'_'+row['ancestor']+'_'+row['intermediate']+'_'+row['leaf']+'_'+row['combined_editing_level_method']+'_strong'+str(row['strong_levels_lower_limit'])+'_'+str(row['strong_levels_upper_limit'])+'_weak'+str(row['weak_levels_lower_limit'])+'_'+str(row['weak_levels_upper_limit'])+'_adapR'+str(row['adaptive_rate'])
                print(fig_name)
                expected_mut,excess_mut,adaptive_mut = [],[],[]
                strong_hp_nonsyn_sites = row['strong_nonsyn_sites']-row['adaptive_sites']
                calc_expected_mutations_distribution(expected_mut,excess_mut,adaptive_mut,
                                                     row['intermediate_syn_a'],row['intermediate_non_syn_a'],row['syn_ag_mut_in_leaf'],row['non_syn_ag_mut_in_leaf'],
                                                     row['weak_syn_sites'],row['weak_nonsyn_sites'],row['strong_syn_sites'],
                                                     strong_hp_nonsyn_sites,row['adaptive_sites'],1,1000000,None,True,False)
                expected_non_syn_eg_mut_dist = [sum(i) for i in zip(expected_mut,excess_mut,adaptive_mut)]
                
                actual_mutations=row['observed_strong_nonsyn_eg_mutations']
                new_p_val = sum([1 for i in expected_non_syn_eg_mut_dist if i<row['observed_strong_nonsyn_eg_mutations']])/float(n)
                row['std']=np.std(expected_non_syn_eg_mut_dist)
                row['p']=new_p_val
                if row['p']==0.0:
                    text='G mut: '+str(int(actual_mutations))+'\nP<1e-6 *'
                else:
                    text='G mut: '+str(int(actual_mutations))+'\nP='+str(pval_str(new_p_val))
                minimal_bin=min([actual_mutations]+expected_non_syn_eg_mut_dist)-5
                maximal_bin=max([actual_mutations]+expected_non_syn_eg_mut_dist)+5
                bins=np.arange(minimal_bin,maximal_bin)
                fig, ax = plt.subplots(figsize=(10,10))
                sns.distplot(expected_non_syn_eg_mut_dist,kde=False,norm_hist=False,bins=bins,color='red')
                # sns.displot(expected_non_syn_eg_mut_dist,kde=False,bins=bins,color='red')
                plt.axvline(x=actual_mutations, color='black', linestyle='--')
                plt.text(0.75,0.05,text,horizontalalignment='left',fontsize=25, transform=ax.transAxes)      
                
                plt.xlim(min([actual_mutations-15,min(expected_non_syn_eg_mut_dist)]),max(max(expected_non_syn_eg_mut_dist),actual_mutations)+10)
                plt.tick_params(axis='both', which='major', labelsize=25)  
                plt.gca().spines['top'].set_visible(False)
                plt.gca().spines['right'].set_visible(False)
                plt.savefig(outpath+fig_name+'.jpg')
                plt.close()
                    
            return row
                    
            
        outpath=path+'expected_mutations_distributions/'
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        results_df=df.apply(lambda row: update_table_and_plot_values(row,name,outpath,n,confidence_level,recalc), axis=1)
        return results_df
    

    #The function starts here!!!!
    results_dfs_dict = {}
    for model_type in models:
            
        try:
            writer = pd.ExcelWriter(path+model_type+'.xlsx', engine='xlsxwriter')
            for iden in iden_list:
                for r in iden_range_list:        
                    name = model_type+'_iden'+str(iden)+'_in_range'+str(r)
                    models_dict = convert_dfs_csv_to_xlsx(path+'*iden'+str(iden)+'_in_range'+str(r)+'_'+model_type)
                    models_dict = convert_models_dfs_back_to_series(models_dict)
                    df=pd.concat([s for s in models_dict.values()], axis=1, sort=False) 
                    df = df.transpose()
                    if recalc:    
                        results_df = restor_and_plot_distributions(path,df,name)
                    else:
                        results_df = df
                    plot_p_values_matrix(path,name,results_df)
                    results_df.to_excel(writer, sheet_name = name)
                    results_dfs_dict.update({'iden'+str(iden)+'_range'+str(r):results_df})
            writer.save()  

        except Exception as e:
            print(e)
            
    return results_dfs_dict
            
        
def calc_dn_ds():
    lst = ['A','T','C','G']
    dn=0
    ds=0
    for c in [''.join(x) for x in it.product(lst, lst, lst)]:
        print('\n')
        print(c)
        if c in ['TAG','TAA','TGA']:
            pass
        else:
            for i,n in enumerate(c):
                if n=="A":
                    new_c = c[:i] + "G" + c[i+1:]
                    print(new_c)
                    if str(Seq(c).translate())==str(Seq(new_c).translate()):
                        ds+=1
                    else:
                        dn+=1
    print(str(ds/(dn+ds)))


def create_tree_relations_dicts(tree_str, is_str=True):
    
    if is_str:
        tree = Phylo.read(StringIO(tree_str), "newick")
    else:
        tree = tree_str
    root_name=tree.root.name
    
    leaves_to_ancestors_dict={}
    leaves=tree.get_terminals()
    for l in leaves:
        ancestors = (root_name,)
        for n in tree.get_path(l.name):
            if n.name != l.name:
                ancestors+=(n.name,)
        leaves_to_ancestors_dict.update({l.name:ancestors})
    
    ancestors_to_leaves_dict={}
    ancestors = tree.get_nonterminals()
    for a in ancestors:
        ancestors_to_leaves_dict.update({a.name:[l.name for l in a.get_terminals()]})
        
    ancestors_to_downstream_ancestors_dict={}
    ancestors = tree.get_nonterminals()
    for a in ancestors:
        ancestors_to_downstream_ancestors_dict.update({a.name:[i.name for i in a.get_nonterminals() if i.name!=a.name]})
    
    ancestors_to_upstream_ancestors_dict={}
    ancestors = tree.get_nonterminals()
    for a in ancestors:
        if a.name==root_name:
            ancestors_to_upstream_ancestors_dict.update({a.name:[]})
        else:
            upstream_ancestors = (root_name,)
            for n in tree.get_path(a.name):
                if n.name!=a.name:
                    upstream_ancestors+=(n.name,)
                ancestors_to_upstream_ancestors_dict.update({a.name:upstream_ancestors})
    
    return tree, leaves_to_ancestors_dict, ancestors_to_leaves_dict, ancestors_to_downstream_ancestors_dict, ancestors_to_upstream_ancestors_dict
                
   
def plt_excess_of_unmutated_sites():
    
    plt.rcdefaults()
    fig, ax = plt.subplots()
    
    # Example data
    tests = ('D to sep', 'D to squ', 'D to bob', 'D to lin', 'S to sep', 'S to squ', 'B to bob', 'B to lin')
    y_pos = np.arange(len(tests))
    performance = [2.32619168,2.217906721,1.909488396,3.744631838,7.626306551,1.757157595,0.471189566,-0.164851777]
    ax.barh(y_pos, performance, align='center')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(tests)
    ax.set_xlim(-2,8)
    ax.invert_yaxis()  # labels read top-to-bottom
    ax.set_xlabel('Excess of unmutated sites\n(expected mut. - actual mut.)')
    plt.axvline(x=0, color='black',linestyle='--')
    plt.show()
    plt.close()
       

def plot_conserved_sites_substitutions(path,rates):
    bars = ()
    y_pos =np.arange(len(bars))
    plt.bar(y_pos, rates, color=['grey','black','red','blue'])
    plt.xticks(y_pos, bars)
    plt.savefig(path+'Gsubs_in_consereved_sites.jpg')
    plt.close()
     

def get_uniqe_proteins_rows(df,animals=['apl','nau','oct','bim','sep','squ']):
    
    def same_protein_name(row,animals):
        proteins = [row[a+'_protein'].split('_')[0] for a in animals]
        if len(set(proteins))==1:
            return True
        else:
            return False
    
    df['same_protein']=df.apply(lambda row: same_protein_name(row,animals), axis=1)
    filtered_df = df[df['same_protein']].copy()
    return filtered_df


def different_aa(row, ancestors_iter):
    if len(set([row[a+'_original_aa'] for a in ancestors_iter]))==1:
        return False
    else:
        return True
    

def create_sites_count_by_animals_table(edited_rows,animals=['oct','bim','sep','squ','bob','lin']):
    
        all_combinations = []
        for r in range(len(animals)+1):
            all_combinations+=list(it.combinations(animals,r))
        
        columns = ['animals','sites_count_for_animals']
        data = []
        for c in all_combinations:
            if len(c):
                rows_for_c = edited_rows
                for a in c:
                    rows_for_c = rows_for_c[rows_for_c[a+'_edited']==1]
                rows_for_c = rows_for_c[rows_for_c['edited_animals']==len(c)]
                data.append(('-'.join(list(c)), len(rows_for_c)))
        
        return pd.DataFrame(data = data, columns = columns)
    
    
def collect_editing_sites(mat, edited_leaves, non_edited_leaves, editing_level_method='average'):
    """
    edited_leaves - a list of tuples. each tuple is a group of animals in which all should be edited
    returns a data frame in which each row comply to at least one of the tuples in edited_leaves
    """
    
    def calc_combined_editing_levels_for_method(row,editing_level_method,edited_animals):
        assert (editing_level_method in ['average','minimal','maximal']),str(editing_level_method)+ "editing level type is not defined, please pass average/minimal/maximal"
        if editing_level_method=='average':
            return float(sum([row[a+'_editing_level'] for a in edited_animals]))/float(sum([row[a+'_edited'] for a in edited_animals]))
        elif editing_level_method=='minimal':
            return min([row[a+'_editing_level'] for a in edited_animals if row[a+'_edited']])
        elif editing_level_method=='maximal':
            return max([row[a+'_editing_level'] for a in edited_animals])


    def in_unification_of_intersections(row, edited_leaves,non_edited_leaves):
        in_unified_intersections = 0
        for group in edited_leaves:
            if len(group)==sum([row[a+'_edited'] for a in group]):
                in_unified_intersections=1
                break
        for leaf in non_edited_leaves:
            if row[leaf+'_edited']==1:
                in_unified_intersections=0
                break
        return in_unified_intersections
    
    edited_animals = list(set([a for group in edited_leaves for a in group]))    
    mat['in_unified_intersections'] = mat.apply(lambda row: in_unification_of_intersections(row, edited_leaves, non_edited_leaves), axis=1)
    mat = mat[mat['in_unified_intersections']==1]  
    editing_level_method=editing_level_method
    edited = mat.drop(labels='in_unified_intersections',axis=1)
    editing_level_method = editing_level_method
    edited['combined_editing_level'] = edited.apply(lambda row: calc_combined_editing_levels_for_method(row, editing_level_method, edited_animals), axis=1)
    return edited
                           
                    
def plot_subs_rates(data,name,path):
    fig, ax = plt.subplots(figsize=(15, 12))
    x = np.arange(len(data))
    
    bar_width = 0.12
    gap = 0
    
    pvals = []
    for d in data:
        pvals += [d[1][4],d[2][4]]  
    pvals_str_list=[]
    for p in pvals:
        if p>0.05:
            pvals_str_list.append('ns')
        elif p<0.0001:
            pvals_str_list.append('****')
        elif p<0.001:
            pvals_str_list.append('***')
        elif p<0.01:
            pvals_str_list.append('**')
        elif p<0.05:
            pvals_str_list.append('*')
    
    x_high = [x[0]-1.5*bar_width-0.5*gap,x[0]+1.5*bar_width+0.5*gap,x[1]-1.5*bar_width-0.5*gap,x[1]+1.5*bar_width+0.5*gap]
    high_labels = ['ancestral A','ancestral G','ancestral A','ancestral G']
    xs = list(x_high)+list(x)
    editing_levels = ['\nWeak sites','\nStrong sites']
    labels = high_labels+editing_levels
    ancestralA_syn = [d[1][0] for d in data]
    ancestralA_syn_err = [d[1][1] for d in data]
    ancestralA_nonsyn = [d[1][2] for d in data]
    ancestralA_nonsyn_err = [d[1][3] for d in data]
    ancestralGsyn = [d[2][0] for d in data]
    ancestralGsyn_err = [d[2][1] for d in data]
    ancestralG_nonsyn = [d[2][2] for d in data]
    ancestralG_nonsyn_err = [d[2][3] for d in data]
    
    b1 = ax.bar(x-2*bar_width+gap, ancestralA_syn, yerr=ancestralA_syn_err, width=bar_width, color='grey', label='syn ancestralA', error_kw=dict(capsize=5))
    b2 = ax.bar(x-bar_width+gap, ancestralA_nonsyn, yerr=ancestralA_nonsyn_err, width=bar_width, color='green', label='non-syn ancestralA', error_kw=dict(capsize=5))
    b3 = ax.bar(x+bar_width-gap, ancestralGsyn, yerr=ancestralGsyn_err, width=bar_width,color='grey', label='syn ancestralG', error_kw=dict(capsize=5))
    b4 = ax.bar(x+2*bar_width-gap, ancestralG_nonsyn, yerr=ancestralG_nonsyn_err, width=bar_width, color='green',label='non-syn ancestralG', error_kw=dict(capsize=5))
    
    
    max_y = max(ancestralA_syn+ancestralA_nonsyn+ancestralGsyn+ancestralG_nonsyn)
    ticks = [-0.3,0.05,0.7,1.05]
    print(pvals_str_list)
    for i,t in enumerate(ticks):   
        plt.plot([t,t+0.25],[max_y*1.25,max_y*1.25],color = 'black')
        ps=pvals_str_list[i]
        if ps=='****':
            align = 0.08
        if ps=='***': 
            align = 0.10
        if ps=='**' or ps=='ns': 
            align = 0.10
        if ps=='*':
            align = 0.11
        plt.text(t+align,max_y*1.27,ps,size=25)
         
    ax.set_xticks(xs)
    ax.set_xticklabels(labels)
    plt.tick_params(axis='both', which='major', labelsize=25)
    ax.set_ylabel('Fraction of sites with G substitutions', labelpad=15,fontsize=34)
    ax.set_ylim([0,max_y*1.35])
    plt.savefig(path+name+'.jpg')
    plt.close()


def calc_substitutions_per_ancestral_nucl(sites_df,levels_ranges,intermediate,animals_to_check_subs,count_subs_multiple_times=True):
    
    rates = []
    if count_subs_multiple_times:
        sites_df['Gsubs'] = sites_df.apply(lambda row: sum([1 for a in animals_to_check_subs if row[a+'_nuc']=="G"]), axis=1)
    else:
        sites_df['Gsubs'] = sites_df.apply(lambda row: 1 if any([row[a+'_nuc']=="G" for a in animals_to_check_subs]) else 0, axis=1)
         
    data = []
    columns = ['editing levels','ancesA_syn','ancesA_syn_subs','ancesA_nonsyn','ancesA_nonsyn_subs','ancesG_syn','ancesG_syn_subs','ancesG_nonsyn','ancesG_nonsyn_subs']
    for el in levels_ranges:
        sites_in_el = sites_df[np.logical_and(sites_df['combined_editing_level']>el[0],sites_df['combined_editing_level']<=el[1])]
        ancestralA = sites_in_el[sites_in_el['N1_nuc']=="A"]
        ancestralG = sites_in_el[sites_in_el['N1_nuc']=="G"]
        ancestralA_syn_n=float(len(ancestralA[ancestralA[intermediate+'_AG_recoding']==0]))
        ancestralA_syn_subs_n=float(ancestralA[ancestralA[intermediate+'_AG_recoding']==0]['Gsubs'].sum())
        ancestralA_nonsyn_n=float(len(ancestralA[ancestralA[intermediate+'_AG_recoding']==1]))
        ancestralA_nonsyn_subs_n=float(ancestralA[ancestralA[intermediate+'_AG_recoding']==1]['Gsubs'].sum())
        ancestralG_syn_n=float(len(ancestralG[ancestralG[intermediate+'_AG_recoding']==0]))
        ancestralG_syn_subs_n=float(ancestralG[ancestralG[intermediate+'_AG_recoding']==0]['Gsubs'].sum())
        ancestralG_nonsyn_n=float(len(ancestralG[ancestralG[intermediate+'_AG_recoding']==1]))
        ancestralG_nonsyn_subs_n=float(ancestralG[ancestralG[intermediate+'_AG_recoding']==1]['Gsubs'].sum())
        ancestralA_syn_rate = ancestralA_syn_subs_n/ancestralA_syn_n
        ancestralA_syn_err = np.sqrt(ancestralA_syn_rate*(1-ancestralA_syn_rate)/ancestralA_syn_n)
        ancestralA_nonsyn_rate = ancestralA_nonsyn_subs_n/ancestralA_nonsyn_n
        ancestralA_nonsyn_err = np.sqrt(ancestralA_nonsyn_rate*(1-ancestralA_nonsyn_rate)/ancestralA_nonsyn_n)
        ancestralG_syn_rate = ancestralG_syn_subs_n/ancestralG_syn_n
        ancestralG_syn_err = np.sqrt(ancestralG_syn_rate*(1-ancestralG_syn_rate)/ancestralG_syn_n)
        ancestralG_nonsyn_rate = ancestralG_nonsyn_subs_n/ancestralG_nonsyn_n
        ancestralG_nonsyn_err = np.sqrt(ancestralG_nonsyn_rate*(1-ancestralG_nonsyn_rate)/ancestralG_nonsyn_n)
        pA=stats.fisher_exact([[ancestralA_syn_n,ancestralA_syn_subs_n],[ancestralA_nonsyn_n,ancestralA_nonsyn_subs_n]])[1]
        pG=stats.fisher_exact([[ancestralG_syn_n,ancestralG_syn_subs_n],[ancestralG_nonsyn_n,ancestralG_nonsyn_subs_n]])[1]
        rates.append((el,(ancestralA_syn_rate,ancestralA_syn_err,ancestralA_nonsyn_rate,ancestralA_nonsyn_err,pA),(ancestralG_syn_rate,ancestralG_syn_err,ancestralG_nonsyn_rate,ancestralG_nonsyn_err,pG)))
        data.append((el,ancestralA_syn_n,ancestralA_syn_subs_n,ancestralA_nonsyn_n,ancestralA_nonsyn_subs_n,ancestralG_syn_n,ancestralG_syn_subs_n,ancestralG_nonsyn_n,ancestralG_nonsyn_subs_n))
    
    data_df = pd.DataFrame(data = data, columns=columns)
    return sites_df, rates, data_df


def calc_substitutions_per_editing_type(sites_df,levels_ranges,intermediate,animals_to_check_subs,count_subs_multiple_times=True):

    def calc_editing_type(row,intermediate,ancestors):
        if row[intermediate+'_AG_recoding']==0:
            return 'syn'
        else:
            if row[intermediate+'_AG_target_aa'] in [row[a+'_original_aa'] for a in ancestors]:
                return 'res'
            else:
                return 'div'
    
    rates = []
    if count_subs_multiple_times:
        sites_df['Gsubs'] = sites_df.apply(lambda row: sum([1 for a in animals_to_check_subs if row[a+'_nuc']=="G"]), axis=1)
    else:
        sites_df['Gsubs'] = sites_df.apply(lambda row: 1 if any([row[a+'_nuc']=="G" for a in animals_to_check_subs]) else 0, axis=1)
        
    sites_df['editing_type'] = sites_df.apply(lambda row: calc_editing_type(row,intermediate,['N1']),axis=1)
        
    data = []
    columns = ['editing levels','syn','syn_subs','res','res_subs','res_p','div','div_subs','div_p','res_div_p']
    for el in levels_ranges:
        sites_in_el = sites_df[np.logical_and(sites_df['combined_editing_level']>el[0],sites_df['combined_editing_level']<=el[1])]
        
        syn = sites_in_el[sites_in_el['editing_type']=="syn"]
        res = sites_in_el[sites_in_el['editing_type']=="res"]
        div = sites_in_el[sites_in_el['editing_type']=="div"]
        syn_n=float(len(syn))
        syn_subs_n=float(syn['Gsubs'].sum())
        res_n=float(len(res))
        res_subs_n=float(res['Gsubs'].sum())
        div_n=float(len(div))
        div_subs_n=float(div['Gsubs'].sum())
        syn_rate = syn_subs_n/syn_n
        syn_err = np.sqrt(syn_rate*(1-syn_rate)/syn_n)
        res_rate = res_subs_n/res_n
        res_err = np.sqrt(res_rate*(1-res_rate)/res_n)
        div_rate = div_subs_n/div_n
        div_err = np.sqrt(div_rate*(1-div_rate)/div_n)
        
        p_res=stats.fisher_exact([[syn_n,syn_subs_n],[res_n,res_subs_n]])[1]
        p_div=stats.fisher_exact([[syn_n,syn_subs_n],[div_n,div_subs_n]])[1]
        p_res_div=stats.fisher_exact([[res_n,res_subs_n],[div_n,div_subs_n]])[1]
        rates.append((el,(syn_rate,syn_err,res_rate,res_err,p_res,div_rate,div_err,p_div,p_res_div)))
        data.append((el,syn_n,syn_subs_n,res_n,res_subs_n,p_res,div_n,div_subs_n,p_div,p_res_div))
    
    data_df = pd.DataFrame(data = data, columns=columns)
    return sites_df, rates, data_df


def get_all_groups_of_edited_sites(tree_str,ancestor,animals=['oct','bim','sep','squ','lin','bob'],rooted=True):

    tree = Phylo.read(StringIO(tree_str), "newick")
    leaves_to_ancestors_dict, ancestors_to_leaves_dict, ancestors_to_downstream_ancestors_dict, ancestors_to_upstream_ancestors_dict = create_tree_relations_dicts(tree) 
    edited_leaves = []
    for c in list(it.combinations(animals,2)):
        common_ances = tree.common_ancestor(c).name
        if common_ances==ancestor or common_ances in ancestors_to_upstream_ancestors_dict[ancestor]:
            edited_leaves.append(c)
    return edited_leaves


def get_groups_of_edited_and_unedited_sites(tree_str,leaf,intermediate,ancestor=None,rooted=True):

    
    tree, leaves_to_ancestors_dict, ancestors_to_leaves_dict, ancestors_to_downstream_ancestors_dict, ancestors_to_upstream_ancestors_dict = create_tree_relations_dicts(tree_str) 
    
    groups_for_product=[]
    other_intermediate_leaves=[l for l in ancestors_to_leaves_dict[intermediate] if l!=leaf and all([k not in set.intersection(set(ancestors_to_downstream_ancestors_dict[intermediate]),set(leaves_to_ancestors_dict[leaf])) for k in leaves_to_ancestors_dict[l]])]
    downstream_to_intermediate_ancestors=ancestors_to_downstream_ancestors_dict[intermediate]
    upstream_to_intermediate_ancestor=[ance for ance in ancestors_to_upstream_ancestors_dict[intermediate] if ance not in ancestors_to_upstream_ancestors_dict[ancestor]]
    for ance in downstream_to_intermediate_ancestors:
        if ance not in leaves_to_ancestors_dict[leaf]:
            groups_for_product.append(ancestors_to_leaves_dict[ance])
    outside_intermediate_leafs=[]
    for ance in upstream_to_intermediate_ancestor:
        outside_intermediate_leafs+=[l for l in ancestors_to_leaves_dict[ance] if l not in ancestors_to_leaves_dict[intermediate]]
    outside_intermediate_leafs=list(set(outside_intermediate_leafs))
    if len(outside_intermediate_leafs):
        groups_for_product.append(outside_intermediate_leafs)
    leaves_for_product=[l for group in groups_for_product for l in group]
    edited_leaves=[]
    for l in other_intermediate_leaves:
        for m in leaves_for_product:
            if l!=m and not len(set.intersection(set(leaves_to_ancestors_dict[l]),set(leaves_to_ancestors_dict[m]),set(ancestors_to_downstream_ancestors_dict[intermediate]))):
                edited_leaves.append(tuple(sorted([l,m])))
    edited_leaves = list(set(edited_leaves))
    edited_leaves = [x for x in edited_leaves if x]
    
    non_edited_leaves = []
    if ancestor is not None:
        all_ancestors_leaves=ancestors_to_leaves_dict[ancestor]
        all_root_leaves=ancestors_to_leaves_dict[tree.root.name]
        for l in all_root_leaves:
            if l not in all_ancestors_leaves:
                non_edited_leaves.append(l)
        non_edited_leaves = non_edited_leaves    
    return edited_leaves, non_edited_leaves
    

def determine_editing_for_ancestor(row,edited_leaves,ancestor):
    """
    for an MSA column row data (pd.Series expects df row)
    Determine if ancestor (str) is edited, 
    based on edited_leaves (expected list of tumples) determined by a tree topology from (e.g. from get_groups_of_edited_and_unedited_sites function)
    """
    if row[ancestor+'_nuc']=="A":
        edited = 0
        for g in edited_leaves:
            if all([row[a+'_edited']==1 for a in g]):
                edited=1
                break
        return edited
    else:
        return 0
    
def clacESmutations(msa_df, tree_str, leaf ,intermediate, mutated_leaf_nucl='G', ancestor=None, rooted=True):
    
    edited_leaves,non_edited_leaves = get_groups_of_edited_and_unedited_sites(tree_str,leaf,intermediate,ancestor=ancestor,rooted=rooted)
    msa_df[intermediate+'_edited']=msa_df.apply(lambda row: determine_editing_for_ancestor(row, edited_leaves, ancestor))
    edited = msa_df[msa_df[intermediate+'_edited']==1]
    mutations = edited[edited[leaf+'_nuc']==mutated_leaf_nucl]
    


def reconstruct_trinity_location(row, animals, trinity_dict, complementary_loc_for_antisense=True):

    for a in animals:
        row[a+'_trinity_component']=row[a+'_id'].split('|')[1]
        trinity_df=trinity_dict[a]
        
        trinity_row = trinity_df.loc[row[a+'_trinity_component'],:]
        strand=trinity_row['strand']
        if row[a+'_coding_location']!='gap':
            
            if strand=="+":
                row[a+'_trinity_locatin_base1']=int(row[a+'_coding_location'])+int(trinity_row['orfs_start'])
            elif strand=="-":
                anti_sense_location_base_1 = int(trinity_row['orfs_end'])-int(row[a+'_coding_location'])-1  
                if complementary_loc_for_antisense:
                    row[a+'_trinity_locatin_base1']=len(trinity_row['sequence'])-anti_sense_location_base_1
                else:
                    row[a+'_trinity_locatin_base1']=anti_sense_location_base_1     
        else:
            row[a+'_trinity_locatin_base1']='-'
    return row
    

def reorganize_sites_data(path,animals,ancestors,tree_strs,trinity_dict):

    df = pd.read_csv(path,sep='\t',index_col=False)
    if 'nau' in animals:
        df = df[~np.logical_and(df['nau_edited']==1, df['edited_animals']==1)]
    df = df.apply(lambda row: reconstruct_trinity_location(row,animals,trinity_dict),axis=1)
    for a in animals:
        df[a+'_aa_swap'] = df.apply(lambda row: row[a+'_original_aa']+row[a+'_AG_target_aa'] if row[a+'_editing_level']>0 else row[a+'_original_aa']+'-', axis=1)
    
    for a in ancestors:
        edited_leaves = get_all_groups_of_edited_sites(tree_str,a)
        df[a+'_edited'] = df.apply(lambda row: determine_editing_for_ancestor(row,edited_leaves,a), axis=1)
        df[a+'_aa_swap'] = df.apply(lambda row: row[a+'_original_aa']+row[a+'_AG_target_aa'] if row[a+'_edited']==1 else row[a+'_original_aa']+'-', axis=1)
        
    df['alignment_position_base1'] = df.apply(lambda row: row['alignment_pos_base_0']+1, axis=1)
    columns = ['super_orthologs_id','alignment_position_base1']
    columns += [a+'_trinity_component' for a in animals]
    columns += [a+'_trinity_locatin_base1' for a in animals]
    columns += [a+'_nuc' for a in animals+ancestors]
    columns += [a+'_editing_level' for a in animals if a not in ['apl','nau']]
    columns += [a+'_aa_swap' for a in animals if a not in ['apl','nau']]
    columns += [a+'_aa_swap' for a in ancestors]
    df = df[columns].copy()
    df.rename(columns={'super_orthologs_id':'orthologs_id','alignment_pos_base_1':'msa_location_base1'}, inplace=True)
    for a in animals+ancestors:
        for c in list(df.columns):
            if a+"_" in c:
                df.rename(columns={c:c.replace(a,animals_names_dict[a])}, inplace=True)
    df.to_excel(path+'_relevant_fields.xlsx',index=False)
    df.to_csv(path+'_relevant_fields.txt',index=False,sep='\t')
    return df


def plot_substitutions_by_type(path, data, types=['synonymous','restorative','diversifying']):
    
    x=np.arange(len(types))
    bar_width=0.5
    fig, ax = plt.subplots(figsize=(15, 12))
    
    syn=data[0]
    syn_err=data[1]
    res=data[2]
    res_err=data[3]
    res_p=data[4] 
    div=data[5] 
    div_err=data[6] 
    div_p=data[7] 
    div_res_p=data[8]
    
    bars = ax.bar(x, [syn,res,div], yerr=[syn_err,res_err,div_err], width=bar_width, color=['grey','deepskyblue','red'], label=types, error_kw=dict(capsize=5))
    ax.set_xticks(x)
    ax.set_ylim([0,0.4])
    ax.set_xticklabels(types)
    ax.set_yticks([0,0.1,0.2,0.3,0.4])
    
    
    txt_size=30
    ax.plot([0,1],[0.3,0.3],color = 'black')
    ax.plot([0,0],[0.3,0.29],color = 'black')    
    ax.plot([1,1],[0.3,0.29],color = 'black')    
    ax.text(0.3,0.31,pval_str(res_p),size=txt_size)
    
    ax.plot([0,2],[0.35,0.35],color = 'black')
    ax.plot([0,0],[0.35,0.34],color = 'black')    
    ax.plot([2,2],[0.35,0.34],color = 'black')    
    ax.text(0.8,0.36,pval_str(div_p),size=txt_size)
    
    ax.plot([1,2],[0.25,0.25],color = 'black')
    ax.plot([1,1],[0.25,0.24],color = 'black')    
    ax.plot([2,2],[0.25,0.24],color = 'black')    
    ax.text(1.3,0.26,pval_str(div_res_p),size=txt_size)
    
    plt.tick_params(axis='x', which='both', length=0, labelsize=30)
    plt.tick_params(axis='y', which='both', labelsize=30)
    ax.set_ylabel('Fraction of sites with G substitutions', labelpad=15,fontsize=35)
    plt.savefig(path+'subs_from_c.jpg')
    plt.close()
       

def plot_nonsyn_to_syn_depletions(animals, rates, errs, path):
    
    labels = animals
    syn_rates = [r[0] for r in rates]
    non_syn_rates = [r[1] for r in rates]
    syn_errs =  [e[0] for e in errs]
    nonsyn_errs = [e[1] for e in errs]
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width/2, syn_rates, width, label='Syn',color='grey',yerr=syn_errs,error_kw=dict(capsize=0))
    rects2 = ax.bar(x + width/2, non_syn_rates, width, label='non-Syn',color='red',yerr=nonsyn_errs,error_kw=dict(capsize=0))
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Editing Events Rate')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend(loc=(0.7,0.6))
    max_y = max([v for r in rates for v in r]) 
    ax.set_ylim([0,max_y*1.3])
    tickes = ax.get_xticks()       
    for i,t in enumerate(tickes):
        plt.plot([t-0.3,t+0.2],[max_y*1.1,max_y*1.1], color="black")
        val = rates[i][1]/rates[i][0]
        val = round_siginficant(val)
        plt.text(t-0.2,max_y*1.12, str(val), color="black")
    plt.savefig(path+'events_rates.jpg')
    plt.close()


def calculates_odds_ratios_z_score(a1,a2,b1,b2,c1,c2,d1,d2):
    a_ratio = a1/float(a2)
    b_ratio = b1/float(b2)
    c_ratio = c1/float(c2)
    d_ratio = d1/float(d2)
    
    ab_odds_ratio = a_ratio/b_ratio
    cd_odds_ratio = c_ratio/d_ratio
    z = (log(ab_odds_ratio)-log(cd_odds_ratio))/sqrt(1/a1+1/a2+1/b1+1/b2+1/c1+1/c2+1/d1+1/d2)
    p = stats.norm.sf(abs(z))
    return z,p


def mutations_rates_tests(path):
    dfs_list = []
    for f in glob.glob(path+'mutations_count*'):
        df = pd.read_csv(f,sep='\t')
        dfs_list.append(df)
    df_all = pd.concat(dfs_list,sort=False)
    df_all['p_val'] = df_all.apply(lambda row: stats.fisher_exact([[row['edited'],row['edited_mutations']],[row['unedited'],row['unedited_mutations']]])[1],axis=1)
    df_all.to_excel(path+'mutations_tests_results.xlsx',index=False)
   


if __name__=='__main__':
    

    
# =============================================================================
#     animals=['oct','bim','sep','squ','bob','lin']
#     conserved_groups = ['bim_oct','sep_squ','bob_lin','bob_lin_sep_squ','bim_bob_lin_oct_sep_squ']
#     path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/all8/hpm_cephalopods_ancestry/'    
#     rates_df, el_df = collect_results_for_hpm_and_plot_probs(path,animals,conserved_groups)      
#     
#     animals=['oct','bim','sep','squ','bob','lin']
#     conserved_groups = ['bim_oct','lin_sep','bob_lin_sep','bob_lin_sep_squ','bim_bob_lin_oct_sep_squ']
#     path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/NCBI/all8/hpm_cephalopods_ancestry_aa/'    
#     rates_df, el_df = collect_results_for_hpm_and_plot_probs(path,animals,conserved_groups)  
# 
#     animals=['oct','bim','sep','squ','bob','lin']
#     conserved_groups = ['bim_oct','lin_squ','bob_lin_squ','bob_lin_sep_squ','bim_bob_lin_oct_sep_squ']
#     path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/oleg/all8/'
#     rates_df, el_df = collect_results_for_hpm_and_plot_probs(path,animals,conserved_groups)  
# =============================================================================
    
    path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/oleg/coleoids_adaptive/'
    results_dfs_dict = collect_results_for_general_model_and_plot_probs(path, recalc=True)
    
    # path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/NCBI/coleoids/new_res/'
    # results_dfs_dict = collect_results_for_general_model_and_plot_probs(path, recalc=True)

    # path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/coleoids/new_res_adaptive/'
    # results_dfs_dict = collect_results_for_general_model_and_plot_probs(path, recalc=True)


    # trinity_files_path='D:/RNA_Editing_large_files_Backup_20201205/transcriptomes_fix/trinity_transcriptomes/our_fix/native_coding_mrna/'
    # animals=['apl','nau','oct', 'bim', 'sep', 'squ', 'lin', 'bob']
    # print('Reading all transcripts from fasta files')
    # trinity_dict = {}
    # for a in animals:
    #     trinity_df = read_fasta_to_df(trinity_files_path+'native_coding_mrna_orfs_'+a+'.fa')
    #     trinity_df.set_index('id', inplace=True)
    #     trinity_dict.update({a:trinity_df})
        
    
    # print('Reading sites file')
    # sites_path = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/coleoids/edited_rows_coleoids'
    # sites_df = pd.read_csv(sites_path,sep='\t',index_col=False)
    # conserved_sets=[('oct',),('bim',),('sep',),('squ',),('bob',),('lin',),
    #                 ('oct','bim'),('sep','squ','bob','lin'),
    #                 ('oct','bim','sep','squ','bob','lin')]
    # outpath= '/'.join(sites_path.split('/')[:-1])+'/codons_dists/'
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath)
    
    # codons_dists_res={}
    # for s in conserved_sets:
    #     print('Calculating codons distributions for '+str(s))
    #     codons_dists_res.update({s:codons_dist_dict(sites_df,trinity_dict,s,aa_instead=False)})

    # dists={}
    # for k,v in codons_dists_res.items():
    #     print('Creating full table for '+str(k))
    #     title='_'.join(k)
    #     mm_df=mismatches_dist(outpath, v, title)
    #     mm_df['editing']=mm_df.apply(lambda row: row['codon']+';'+str(row['nucl']),axis=1)
    #     mm_df.set_index('editing', inplace=True)
    #     dists.update({k:mm_df})
    
    # animals=['oct', 'bim', 'sep', 'squ', 'lin', 'bob']
    # print('merging tables for '+str(animals))
    # df_final=dists[(animals[0],)]
    # df_final.rename(columns={animals[0]+'_events':animals_names_dict[animals[0]]}, inplace=True)
    # for a in animals[1:]:
    #     df_final[animals_names_dict[a]]=dists[(a,)][a+'_events']
    # df_final.to_excel(outpath+'edited_codons.xlsx', index=False)
    # title='edited_codons_dist'
    # plot_mm_dist_for_multiple_animals(df_final,outpath,title,[animals_names_dict[a] for a in animals])
    
    
    
    
        
        
# =============================================================================
#     neural_path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/dnds_analysis/dnds_neural'
#     non_neural_path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/dnds_analysis/dnds_non_neural'
#     df1 = pd.read_csv(neural_path,sep='\t',index_col=0)
#     df2 = pd.read_csv(non_neural_path,sep='\t',index_col=0)
# =============================================================================
    
    
# =============================================================================
#     s = "O.vult19956t681198t36527t1651609\nO.bimt12868t687196t23751t1664678\nS.ofit21942t650796t41986t1594042\nD.peat13816t684107t25836t1637114\nEupt8716t668266t16140t1651819\nS.lint16327t675324t29367t1621529\n"
#     rows = s.split('\n')
#     rows=rows[:-1]
#     data=[]
#     for r in rows:
#         row_data=r.split('t')
#         row_data = [row_data[0]]+[int(n) for n in row_data[1:]]
#         data.append(row_data)
#     df=pd.DataFrame(data=data,columns=['animal','syn_sites','syn_a','nonsyn_sites','nonsyn_a'])
#     animals=list(df.animal.values)
#     df.set_index('animal', inplace=True)
#     
#     rates=[]
#     errs=[]
#     for a in animals:
#         
#         syn_r = float(df.loc[a,'syn_sites'])/float(df.loc[a,'syn_a'])
#         nonsyn_r = float(df.loc[a,'nonsyn_sites'])/float(df.loc[a,'nonsyn_a'])
#         syn_err = np.sqrt(syn_r*(1-syn_r)/df.loc[a,'syn_a'])
#         nonsyn_err = np.sqrt(nonsyn_r*(1-nonsyn_r)/df.loc[a,'nonsyn_a'])
#         rates.append((syn_r,nonsyn_r))
#         errs.append((syn_err,nonsyn_err))
#         
#     outpath = 'E:/RNA_editing_Large_files/Phylogeny/results/coleoids/our_model/rooted/'
#     plot_nonsyn_to_syn_depletions(animals, rates, errs, outpath)
# =============================================================================
    
    
# =============================================================================
#     trinity_files_path='C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/transcriptomes_fix/trinity_transcriptomes/our_fix/new_native_transcriptomes/'
#     animals=['apl','nau','oct', 'bim', 'sep', 'squ', 'lin', 'bob']
#     print('Reading all transcripts from fasta files')
#     trinity_dict = {}
#     for a in animals:
#         trinity_df = read_trinity_mrna_files(trinity_files_path+'new_native_orfs_'+a+'.fa')
#         trinity_df.set_index('id', inplace=True)
#         trinity_dict.update({a:trinity_df})
#     
#     path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/coleoids/edited_rows_coleoids'
#     tree_str = trees['coleoids_rooted']
#     print("Processing table for \n"+tree_str)
#     animals=['oct', 'bim', 'sep', 'squ', 'lin', 'bob']
#     ancestors=['C','O','D','B','S']
#     df_coleoids=reorganize_sites_data(path,animals,ancestors,tree_str,trinity_dict)
#     
#     path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/all8/edited_rows_all8'
#     tree_str = trees['all8_rooted']
#     print("Processing table for \n"+tree_str)
#     animals=['apl','nau','oct', 'bim', 'sep', 'squ', 'lin', 'bob']
#     ancestors=['N0','N1','C','O','D','B','S']
#     df_all8=reorganize_sites_data(path,animals,ancestors,tree_str,trinity_dict)
# =============================================================================
    
# =============================================================================
# #    path = 'E:/RNA_editing_Large_files/Phylogeny/results/coleoids/edited_rows_coleoids'
#     path = 'E:/RNA_editing_Large_files/Phylogeny/results/all8/edited_rows_all8'
#     animals = ['oct','bim','sep','squ','bob','lin']
#     ancestors = ['C','O','D','S','B']
#     columns = ['super_orthologs_id','alignment_pos_base_0']
#     columns += [a+'_nuc' for a in animals+ancestors]
#     columns += [a+'_edited' for a in animals]
#     columns += [a+'_editing_level' for a in animals]
#     columns += [a+'_original_aa' for a in animals+ancestors]
#     columns += [a+'_AG_recoding' for a in animals+ancestors]
#     columns += [a+'_AG_target_aa' for a in animals+ancestors]
#     df = pd.read_csv(path,sep='\t',index_col=False)
#     df = df[columns].copy()
#     df.to_excel(path+'.xlsx',index=False)
#     df.to_csv(path+'_relevant_fields.txt',index=False,sep='\t')
# =============================================================================
    
# =============================================================================
#     path = 'E:/RNA_editing_Large_files/Phylogeny/results/all8/edited_rows'
#     animals = ['apl','nau','oct','bim','sep','squ','bob','lin']
#     ancestors = ['N0','N1','C','O','D','S','B']
#     columns = ['super_orthologs_id','alignment_pos_base_0']
#     columns += [a+'_nuc' for a in animals+ancestors]
#     columns += [a+'_edited' for a in animals]
#     columns += [a+'_editing_level' for a in animals]
#     columns += [a+'_original_aa' for a in animals+ancestors]
#     columns += [a+'_AG_recoding' for a in animals+ancestors]
#     columns += [a+'_AG_target_aa' for a in animals+ancestors]
#     
#     df = pd.read_csv(path,sep='\t',index_col=False)
#     df = df[columns].copy()
#     df.to_excel(path+'.xlsx',index=False)
#     
#     print(str(len(df[df['N1_nuc']=='A'])) + ' N1 ancestralA')
#     print(str(len(df[df['N1_nuc']=='G'])) + ' N1 ancestralG')
#     print(str(len(df[df['N1_nuc']=='C'])) + ' N1 ancestralC')
#     print(str(len(df[df['N1_nuc']=='T'])) + ' N1 ancestralT')
# =============================================================================
    
    
# =============================================================================
#     # params_groups = [
#     #                  ('C',([0,0.1],[0.1,1],[0,1]),[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],['sep','squ','bob','lin','oct','bim'],'all8_rooted',
#     #                  'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/non_neural/4fold/edited_rows_from_4fold_non_neural_subset','non_neural'),
#     #                  ('C',([0,0.1],[0.1,1],[0,1]),[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],['sep','squ','bob','lin','oct','bim'],'all8_rooted',
#     #                  'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/neural/4fold/edited_rows_from_4fold_neural_subset','neural')
#                      
#     #                  ]
# 
#     params_groups = [
#                      ('C',([0,0.1],[0.1,1],[0,1]),[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],['sep','squ','bob','lin','oct','bim'],'all8_rooted_ncbi',
#                       'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/all8/edited_rows','raxml'),
#         
#                      ('C',([0,0.1],[0.1,1],[0,1]),[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],['sep','squ','bob','lin','oct','bim'],'all8_rooted_oleg',
#                      'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/NCBI/all8/edited_rows','oleg'),
#                      
#                      ('C',([0,0.1],[0.1,1],[0,1]),[('sep','oct'),('sep','bim'),('squ','oct'),('squ','bim'),('bob','oct'),('bob','bim'),('lin','oct'),('lin','bim')],['sep','squ','bob','lin','oct','bim'],'all8_rooted_ncbi',
#                       'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/oleg/all8/edited_rows','ncbi')
#                      ]
# 
#     
#     results = []
#     conserved_sites_subs_dict = {}
#     for param in params_groups:
#         
#         intermediate = param[0] 
#         levels_ranges = param[1]
#         edited_leaves = param[2]
#         animals_to_check_subs=param[3]
#         tree=param[4]
#         path=param[5]
#         params_name=param[6]
#         df = pd.read_csv(path,sep='\t',index_col=False)
#         non_edited_leaves = []
#         df_for_intermediate = df[df[intermediate+'_nuc']=="A"]
#         relevant_es = collect_editing_sites(df_for_intermediate, edited_leaves, non_edited_leaves, editing_level_method='average')
#         print(str(len(relevant_es[relevant_es['N1_nuc']=='A'])) + ' N1 ancestralA')
#         print(str(len(relevant_es[relevant_es['N1_nuc']=='G'])) + ' N1 ancestralG')
#         print(str(len(relevant_es[relevant_es['N1_nuc']=='C'])) + ' N1 ancestralC')
#         print(str(len(relevant_es[relevant_es['N1_nuc']=='T'])) + ' N1 ancestralT')
#         sites_df, rates,data_df = calc_substitutions_per_ancestral_nucl(relevant_es,levels_ranges,intermediate,animals_to_check_subs,count_subs_multiple_times=False)
#         sites_df, rates, data_df = calc_substitutions_per_editing_type(relevant_es,levels_ranges,intermediate,animals_to_check_subs,count_subs_multiple_times=False)
#         
#         results.append(rates)
#         outpath = '/'.join(path.split('/')[:-1])+'/'
#         name = 'subs_from_'+intermediate
#         data_df.to_excel(path+name+'.xlsx', index=False)
#         plot_substitutions_by_type(outpath,rates[2][1])
#         conserved_sites_subs_dict.update({params_name:data_df})
# =============================================================================
 
    
    
# =============================================================================
#     path = 'E:/RNA_editing_Large_files/Phylogeny/results/no_boblin/unrooted/hpm/edited_at_least_3'
#     outpath = '/'.join(path.split('/')[:-1])+'/'+path.split('/')[-1]+'.xlsx'
#     edited_at_least_3 = pd.read_csv(path, sep='\t', index_col=False)
#     columns_to_write = [col for col in list(edited_at_least_3.columns) if 'id' in col]+['alignment_pos_base_0']
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if '_nuc' in col]
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if '_edited' in col]
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if 'AG_recoding' in col]
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if 'editing_level' in col]
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if 'original_aa' in col]
#     columns_to_write += [col for col in list(edited_at_least_3.columns) if 'AG_target_aa' in col]
#     columns_to_write += ['has_a_genomicG','editing_type']
#     for c in list(edited_at_least_3.columns):
#         if c not in columns_to_write:
#             edited_at_least_3 = edited_at_least_3.drop(c, axis=1)
# #    edited_at_least_3 = edited_at_least_3.loc[columns_to_write]
#     edited_at_least_3.to_excel(outpath,index=False)
#     
#     conserved_4_coleoids = pd.read_csv(path, sep='\t', index_col=False)
#     conserved_4_coleoids_same_prots = get_uniqe_proteins_rows(conserved_4_coleoids)
#     conserved_4_coleoids_same_prots['gene_name']  = conserved_4_coleoids_same_prots.apply(lambda row: row['oct_protein'].split('_')[0], axis =1)
#     conserved_4_coleoidsGsubs = conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['has_a_genomicG']==1]    
#     prots_in_results = list(set(conserved_4_coleoids_same_prots['gene_name']))
#  
#     syn=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']=='syn'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']=='syn'])
#     div=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']=='div'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']=='div'])
#     non_syn=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']!='syn'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']!='syn'])
#     plot_conserved_sites_substitutions(path,[syn,div,non_syn])
# =============================================================================
      

    
# =============================================================================
#     path = 'C:/Users/shosh/OneDrive/Desktop/results/all8/hpm/edited_at_least_3'
#     conserved_4_coleoids = pd.read_csv(path, sep='\t', index_col=False)
#     conserved_4_coleoids_same_prots = get_uniqe_proteins_rows(conserved_4_coleoids)
#     conserved_4_coleoids_same_prots['gene_name']  = conserved_4_coleoids_same_prots.apply(lambda row: row['oct_protein'].split('_')[0], axis =1)
#     conserved_4_coleoidsGsubs = conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['has_a_genomicG']==1]    
#     prots_in_results = list(set(conserved_4_coleoids_same_prots['gene_name']))
# 
#     syn=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']=='syn'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']=='syn'])
#     div=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']=='div'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']=='div'])
#     non_syn=len(conserved_4_coleoidsGsubs[conserved_4_coleoidsGsubs['editing_type']!='syn'])/len(conserved_4_coleoids_same_prots[conserved_4_coleoids_same_prots['editing_type']!='syn'])
#     plot_conserved_sites_substitutions(path,[syn,div,non_syn])
# =============================================================================

# =============================================================================
#     path = 'C:/Users/shosh/OneDrive/Desktop/results/all8/'
#     ances_mismatch_AG = 'N0A_N1G
#     file = 'C:/Users/shosh/OneDrive/Desktop/results/all8/' + ances_mismatch_AG
#     unmatching_ances_AG = pd.read_csv(file, sep='\t', index_col=False)
#     unmatching_ances_AG['different_aa'] = unmatching_ances_AG.apply(lambda row: different_aa(row, ['N0','N1']), axis=1)
#     unmatching_ances_AG = unmatching_ances_AG[~np.logical_or(unmatching_ances_AG['apl_nuc']=='gap',np.logical_or(unmatching_ances_AG['nau_nuc']=='gap',unmatching_ances_AG['gapped_animals']>3))]
#     print('unmatching nucl (no gaps in apl/nau/>3): '+ str(len(unmatching_ances_AG)))
#     unmatching_ances_AG=unmatching_ances_AG[unmatching_ances_AG['different_aa']]
#     print('different aa: '+str(len(unmatching_ances_AG)))
#     print('editing sites: '+str(len(unmatching_ances_AG[unmatching_ances_AG['edited_animals']>0])))
#     print('conserved editing sites: '+str(len(unmatching_ances_AG[unmatching_ances_AG['edited_animals']>1])))
#     aa_columns_AG=unmatching_ances_AG[['apl_bim_bob_lin_nau_oct_sep_squ_consensus_ratio_10_aa_range','super_orthologs_id','alignment_length','alignment_pos_base_0','edited_animals','gapped_animals','A_animals','C_animals','T_animals','G_animals']+[a+'_nuc' for a in ['N0','N1','C','O','D','S','B','apl','nau','oct','bim','sep','squ','bob','lin']]+[a+'_original_aa' for a in ['N0','N1','C','O','D','S','B','apl','nau','oct','bim','sep','squ','bob','lin']]]
#     print('unmatching nucl (no gaps): '+str(len(aa_columns_AG[aa_columns_AG['gapped_animals']==0])))
#     len(aa_columns_AG[aa_columns_AG['apl_nuc']=='C'])
#     print(len(aa_columns_AG[aa_columns_AG['apl_nuc']=='A']))
#     aa_columns_AG_no_aplA = aa_columns_AG[aa_columns_AG['apl_nuc']!='A']
#     
#     path = 'C:/Users/shosh/OneDrive/Desktop/results/all8/'
#     ances_mismatch_GA = 'N0G_N1A'
#     file = 'C:/Users/shosh/OneDrive/Desktop/results/all8/' + ances_mismatch_GA
#     unmatching_ances_GA = pd.read_csv(file, sep='\t', index_col=False)
#     unmatching_ances_GA['different_aa'] = unmatching_ances_GA.apply(lambda row: different_aa(row, ['N0','N1']), axis=1)
#     unmatching_ances_GA = unmatching_ances_GA[~np.logical_or(unmatching_ances_GA['apl_nuc']=='gap',np.logical_or(unmatching_ances_GA['nau_nuc']=='gap',unmatching_ances_GA['gapped_animals']>3))]
#     print('unmatching nucl (no gaps in apl/nau/>3): '+ str(len(unmatching_ances_GA)))
#     unmatching_ances_GA=unmatching_ances_GA[unmatching_ances_GA['different_aa']]
#     print('different aa: '+str(len(unmatching_ances_GA)))
#     print('editing sites: '+str(len(unmatching_ances_GA[unmatching_ances_GA['edited_animals']>0])))
#     print('conserved editing sites: '+str(len(unmatching_ances_GA[unmatching_ances_GA['edited_animals']>1])))
#     aa_columns_GA=unmatching_ances_GA[['apl_bim_bob_lin_nau_oct_sep_squ_consensus_ratio_10_aa_range','super_orthologs_id','alignment_length','alignment_pos_base_0','edited_animals','gapped_animals','A_animals','C_animals','T_animals','G_animals']+[a+'_nuc' for a in ['N0','N1','C','O','D','S','B','apl','nau','oct','bim','sep','squ','bob','lin']]+[a+'_original_aa' for a in ['N0','N1','C','O','D','S','B','apl','nau','oct','bim','sep','squ','bob','lin']]]
#     print('unmatching nucl (no gaps): '+str(len(aa_columns_GA[aa_columns_GA['gapped_animals']==0])))
#     print(len(aa_columns_GA[aa_columns_GA['apl_nuc']=='G']))
#     aa_columns_GA_no_aplG = aa_columns_GA[aa_columns_GA['apl_nuc']!='G']
# =============================================================================
    
    
    
# =============================================================================
# #    path = 'E:/RNA_editing_Large_files/Phylogeny/results/coleoids/our_model/rooted/edited_rows'
# #    animals=['oct','bim','sep','squ','bob','lin']
#     path = 'E:/RNA_editing_Large_files/Phylogeny/results/all8/our_model/rooted/edited_rows'
#     animals=['oct','bim','sep','squ','bob','lin']
#     
#     edited_rows = pd.read_csv(path, sep='\t',index_col=False)
#     df_edited_sites_cnt = create_sites_count_by_animals_table(edited_rows,animals=animals)
#     out_path = '/'.join(path.split('/')[:-1])+'/'
#     df_edited_sites_cnt.to_excel(out_path+'species_sites_count.xlsx',index=False)
# =============================================================================
    
    

# =============================================================================
#     # General transctiprom and editing sites data tables
#     animals = ['apl','nau','oct', 'bim', 'sep', 'squ', 'bob', 'lin']
#     edited_animals = ['oct', 'bim', 'sep', 'squ', 'bob', 'lin']
#     transcripts_path = 'E:/RNA_editing_Large_files/transcriptomes_fix/trinity_transcriptomes/our_fix/new_native_transcriptomes/'
#     sites_path = 'E:/RNA_editing_Large_files/transcriptomes_fix/trinity_transcriptomes/our_fix/filtered_editing_sites/'
#     outpath  = 'C:/Users/shosh/Google_Drive/RNA_Editing/Phylogeny/Tables/'
# 
#     orfs_columns = ('Species','Number of ORFs','Number of unique proteins','Mean ORF length (nt)','Median ORF length (nt)','Total ORF length (nt)')
#     sites_columns = ('Species','AG sites','recoding sites','ORFs containing sites', 'ORFs containing recoding sites')
#     orfs_data = []
#     sites_data = []
#     
#     for a in animals:
#         
#         trinity_df = read_trinity_mrna_files(transcripts_path+'new_native_orfs_'+a+'.fa')    
#         trinity_df['length'] = trinity_df.apply(lambda row: row['orfs_end']-row['orfs_start'],axis=1)
#         orfs_n = len(trinity_df)
#         unique_prot = len(set(list(trinity_df['protein'].values)))
#         mean_l = np.mean(trinity_df['length'].values)
#         median_l = np.median(trinity_df['length'].values)
#         total_l = sum(trinity_df['length'].values)    
#         orfs_data.append((animals_names_dict[a],orfs_n,unique_prot,mean_l,median_l,total_l))
#         
#         if a in edited_animals:
#             sites_df = read_editing_sites_tbl(sites_path+a+'_editing_site_info.txt') 
#             sites_df['recoding']=sites_df.apply(lambda row: 1 if row['AA_before']!=row['AA_after'] else 0, axis=1)
#             sites_n = len(sites_df)
#             sites_recoding_n = len(sites_df[sites_df['recoding']==1])
#             orfs_n = len(set(list(sites_df['id'].values)))
#             orfs_recoding_n = len(set(list(sites_df[sites_df['recoding']==1]['id'].values)))
#             sites_data.append((animals_names_dict[a],sites_n,sites_recoding_n,orfs_n,orfs_recoding_n))
#     
#     sites_data_df = pd.DataFrame(data=sites_data,columns=sites_columns)
#     sites_data_df.to_excel(outpath+'sites_data.xlsx',index=False)
#     orfs_data_df = pd.DataFrame(data=orfs_data,columns=orfs_columns)
#     orfs_data_df.to_excel(outpath+'orfs_data.xlsx',index=False)
# =============================================================================

