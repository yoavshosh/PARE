# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 13:50:41 2020

@author: shosh
"""

import pandas as pd
import numpy as np
import os
import sys
import glob
import pickle
from plots_for_hypotheses_test_results import *

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import lines 
    import matplotlib.transforms as transforms
    import pylab
except ImportError:
    print('Could not Import data-viz libraries')


def calc_average(row,tissue):
    return sum([row[c] for c in list(final_tbl_only_tpm.columns) if tissue in c])/float(len([row[c] for c in list(final_tbl_only_tpm.columns) if tissue in c]))

    
def dnds_odds_ratio_analysis(df1,df2, name1='neural',name2='non_neural'):
    df1['dnds'] = df1.apply(lambda row: float(row['nonsyn_mutated'])*row['syn_nucl']/(row['nonsyn_nucl']*float(row['syn_mutated'])), axis=1)
    df2['dnds'] = df2.apply(lambda row: float(row['nonsyn_mutated'])*row['syn_nucl']/(row['nonsyn_nucl']*float(row['syn_mutated'])), axis=1)
    merged_df = df1.merge(df2,suffixes=('_'+name1,'_'+name2),left_index=True,right_index=True)
    merged_df['Z'] = merged_df.apply(lambda row: calculates_odds_ratios_z_score(row['nonsyn_mutated_neural'],row['nonsyn_nucl_neural'],row['syn_mutated_neural'],row['syn_nucl_neural'],row['nonsyn_mutated_non_neural'],row['nonsyn_nucl_non_neural'],row['syn_mutated_non_neural'],row['syn_nucl_non_neural'])[0], axis=1)
    merged_df['pval'] = merged_df.apply(lambda row: calculates_odds_ratios_z_score(row['nonsyn_mutated_neural'],row['nonsyn_nucl_neural'],row['syn_mutated_neural'],row['syn_nucl_neural'],row['nonsyn_mutated_non_neural'],row['nonsyn_nucl_non_neural'],row['syn_mutated_non_neural'],row['syn_nucl_non_neural'])[1], axis=1)
    return merged_df

def plot_odds_ratios_for_sites_groups(outfile,hpm_rates_dict, animals=['oct','bim','sep','squ','bob','lin']):
    
    data = []
    columns = ('Species','syn_neural','syn_neural_err','div_neural','div_neural_err','syn_non_neural','syn_non_neural_err','div_non_neural','div_non_neural_err',)
    
    p_vals = {}
    odds_ratios = {}
    
    for a in animals:
        a1 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'div_edited']
        a2 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'div']
        b1 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'syn_edited']
        b2 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'syn']
        c1 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'div_edited']
        c2 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'div']
        d1 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'syn_edited']
        d2 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'syn']
        z,p = calculates_odds_ratios_z_score(a1,a2,b1,b2,c1,c2,d1,d2)
        
        syn_neural = float(b1)/b2
        div_neural = float(a1)/a2
        syn_non_neural = float(d1)/d2
        div_non_neural = float(c1)/c2
        syn_neural_err = np.sqrt(syn_neural*(1-syn_neural)/b2)
        div_neural_err = np.sqrt(div_neural*(1-div_neural)/a2)
        syn_non_neural_err = np.sqrt(syn_non_neural*(1-syn_non_neural)/d2)
        div_non_neural_err = np.sqrt(div_non_neural*(1-div_non_neural)/c2)
        
        data.append(([animals_names_dict[a],syn_neural,syn_neural_err,div_neural,div_neural_err,syn_non_neural,syn_non_neural_err,div_non_neural,div_non_neural_err]))
        
        odds_ratios.update({animals_names_dict[a]:(div_neural/syn_neural, div_non_neural/syn_non_neural)})
        p_vals.update({animals_names_dict[a]:pval_str(p)})
        
    df = pd.DataFrame(data=data,columns=columns)
    df.set_index('Species',inplace=True)
    
    max_y = max(list(df['syn_neural'])+list(df['div_neural'])+list(df['syn_non_neural'])+list(df['div_non_neural']))        
    pdax = pd.concat([df.syn_neural.rename('syn-neural'), df.div_neural.rename('div-neural'), df.syn_non_neural.rename('syn-non neural'),df.div_non_neural.rename('div-non neural')], axis=1).plot(kind='bar', yerr=(df['syn_neural_err'].values,df['div_neural_err'].values,df['syn_non_neural_err'].values,df['div_non_neural_err'].values), error_kw=dict(lw=1, capsize=2, capthick=1), color=['grey','red','lightgrey','lightsalmon'])
    ticks = pdax.get_xticks()
    labels = pdax.get_xticklabels()
    plt.xticks(rotation=0)
    plt.title('Diversifying and Synonymous Sites\nin Neural and Non-Neural Genes'+'\n\n', fontsize=17)
    plt.ylabel('Fraction of sites edited', fontsize=17)
    plt.xlabel('Species', fontsize=17)
    plt.tick_params(axis='x', which='both', length=0, labelsize=13)
    plt.tick_params(axis='y', which='both', labelsize=15)
    plt.ylim([0,max_y*1.3])
    plt.gcf().set_size_inches(10,7)
    pdax.legend(ncol=4,fontsize=13,bbox_to_anchor=(1.01, 1.1))
    for i,t in enumerate(ticks):
        plt.plot([t-0.35,t+0.4],[max_y*1.2,max_y*1.2],color = 'black')
        plt.plot([t-0.35,t-0.35],[max_y*1.2,max_y*1.18],color = 'black')
        plt.plot([t+0.4,t+0.4],[max_y*1.2,max_y*1.18],color = 'black')
        ps=p_vals[labels[i].get_text()]
        if ps=='ns':
            align=-0.1
        else:
            align=-0.25
        plt.text(t+align,max_y*1.22,ps,size=10)
        plt.plot([t-0.32,t+0],[max_y*1.1,max_y*1.1],color = 'black')
        plt.plot([t-0.32,t-0.32],[max_y*1.1,max_y*1.08],color = 'black')
        plt.plot([t+0,t+0],[max_y*1.1,max_y*1.08],color = 'black')
        dnds=str(round_siginficant(odds_ratios[labels[i].get_text()][0]))
        if float(dnds)>=1:
            align=-0.25
        else:
            align=-0.29
        plt.text(t+align,max_y*1.12,dnds,size=10)
        plt.plot([t+0.06,t+0.38],[max_y*1.1,max_y*1.1],color = 'black')
        plt.plot([t+0.38,t+0.38],[max_y*1.1,max_y*1.08],color = 'black')
        plt.plot([t+0.06,t+0.06],[max_y*1.1,max_y*1.08],color = 'black')
        dnds=str(round_siginficant(odds_ratios[labels[i].get_text()][1]))
        if float(dnds)>=1:
            align=0.1
        else:
            align=0.08
        plt.text(t+align,max_y*1.12,dnds,size=10)
        
    plt.savefig(outfile)
    plt.close()
        
    
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx 



if __name__=='__main__':
    
    create=True
    plot=True
    
    if create: # FOR CREATING SUBSEETS OF NEURAL/NON-NEURAL GENES BASED ON SALMON RESULTS PER TISSUE
    
        neural_tiss = ['optic_lobe','axial nerve_cord','subesophageal_brain','supraesophageal_brain']
        non_neural_tiss = ['ova','testes','posterior_salivary_gland','skin','suckers','embryos','retina','hepatopancreas_kidney_heart']
        quant_columns = ['Name','Length','EffectiveLength','TPM','NumReads']
        
        # salmon_quant_path = sys.argv[1]
        # sra_metadata_tbl = sys.argv[2]
        # sites_tbl = sys.argv[3]
        # out_dir = os.getcwd()+'/'+sys.argv[4]+'/'
        # non_neural_per_neural=int(sys.argv[5])
            
        salmon_quant_path = 'D:/RNA_Editing_large_files_Backup_20201205/SALMON/quant/'
        sra_metadata_tbl = 'D:/RNA_Editing_large_files_Backup_20201205/SALMON/SraRunTable.txt'
        sites_tbl = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Raxml/coleoids/edited_rows_coleoids'
        out_dir = '/'.join(sra_metadata_tbl.split('/')[:-1])+'/classifications_per_tissue_type/'
        non_neural_per_neural = 1
    
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        # Process Salmon Results for different samples    
        print('Reading sra_metadata_tbl')
        sra_tbl = pd.read_csv(sra_metadata_tbl,sep=',',engine='python')
        sra_tbl.set_index('Run',inplace=True)
        sra_names = []    
        sra_dfs = []
        quant_results=glob.glob('/'.join([salmon_quant_path,'*_quant/quant.sf']))
        print('Constructing SALMON results tbl')
        for f in quant_results:
            names_dict = {}
            f=f.replace('\\','/')
            sra_name = f.split('/')[-2].split('_')[0]
            sra_names.append(sra_name)
            tissue = sra_tbl.loc[sra_name,'Tissue'].replace(' ','_').replace('\\','').replace(',','')
            [names_dict.update({c:c+'_'+sra_name+'_'+tissue}) for c in quant_columns]
            df = pd.read_csv(f,sep='\t',index_col=False)
            df.rename(columns=names_dict,inplace=True)
            df.set_index('Name'+'_'+sra_name+'_'+tissue,inplace=True)
            sra_dfs.append(df)
        
        final_tbl = pd.concat(sra_dfs,axis=1)
        final_tbl_only_tpm = final_tbl[[c for c in final_tbl.columns if "TPM" in c]].copy()
        
        print('Calculating results per tissue type')
        for t in sra_tbl['Tissue'].values:
            tiss = t.replace(' ','_').replace('\\','').replace(',','')
            final_tbl_only_tpm[tiss+'_average']=final_tbl_only_tpm.apply(lambda row: calc_average(row,tiss), axis=1)
        
        final_tbl_only_tpm_averaged  = final_tbl_only_tpm[[c for c in final_tbl_only_tpm if 'average' in c]].copy()
        final_tbl_only_tpm_averaged['neural_average']=final_tbl_only_tpm_averaged.apply(lambda row: 
                                                                                        sum([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in neural_tiss])])/
                                                                                        float(len([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in neural_tiss])])), axis=1)
        
        final_tbl_only_tpm_averaged['non_neural_average']=final_tbl_only_tpm_averaged.apply(lambda row: 
                                                                                        sum([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in non_neural_tiss])])/
                                                                                        float(len([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in non_neural_tiss])])), axis=1)
        
        final_tbl_only_tpm_averaged['neural_xfold'] = final_tbl_only_tpm_averaged.apply(lambda row: row['neural_average']/row['non_neural_average'] if row['non_neural_average']!=0 else np.inf, axis=1)
        final_tbl_only_tpm_averaged.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/SALMON_tmp_results_all_samples_averaged_tpm.xlsx')
        
        
    
        
        # Define and save neural and non-neural subsets for different tresholds
        final_tbl_only_tpm_averaged = final_tbl_only_tpm_averaged.reset_index()
        final_tbl_only_tpm_averaged['index']=final_tbl_only_tpm_averaged.apply(lambda row: 'bim|'+row['index'], axis=1)
        final_tbl_only_tpm_averaged.set_index('index',inplace=True)
        
        
        print('collecting sets from DB for different cutoofs ')
        sites_df = pd.read_csv(sites_tbl,sep='\t', index_col=False)
        final_tbl_only_tpm_averaged_reduced=final_tbl_only_tpm_averaged.loc[list(set(list(sites_df['bim_id'].values))),:].copy()
        
        tpm_rnage=0.1
        for x in [1.5,2,2.5,3,4,5,7,10]:
            
            print('\n\nTesting for '+str(x)+'-fold')
            #collect neural genes and sites based on x-fold expression
            neural_genes_subset = final_tbl_only_tpm_averaged_reduced[np.logical_and(final_tbl_only_tpm_averaged_reduced['neural_xfold']>=x, final_tbl_only_tpm_averaged_reduced['neural_average']>1)].index
            neural_sites_subset_from_sites = sites_df[sites_df['bim_id'].isin(neural_genes_subset)]
            
            #smae for non neural
            non_neural_subset_df = final_tbl_only_tpm_averaged_reduced[~final_tbl_only_tpm_averaged_reduced.index.isin(neural_genes_subset)]
            non_neural_genes_subset = list(non_neural_subset_df.index)
            non_neural_subset_tpms = list(non_neural_subset_df['non_neural_average'].values)
            non_neural_sites_subset_from_sites = sites_df[sites_df['bim_id'].isin(non_neural_genes_subset)]
            
            #collect non-neural subset based on neural genes tpm - for each neural gene draw a non-neural genes with similar tpm in non-neural tissues
            non_neural_per_neural_genes_subset=[]
            diff_lst = []
            for neural_gene in neural_genes_subset:
                neural_tpm = final_tbl_only_tpm_averaged_reduced.loc[neural_gene,'neural_average']
                for nn in range(non_neural_per_neural):
                    nearest_non_neural_tpm_index = find_nearest(non_neural_subset_tpms,neural_tpm)
                    non_neural_tpm=non_neural_subset_tpms[nearest_non_neural_tpm_index]
                    diff=abs(non_neural_tpm-neural_tpm)/float(neural_tpm)
                    diff_lst.append(diff)
                    non_neural_gene=non_neural_genes_subset[nearest_non_neural_tpm_index]
                    if diff>tpm_rnage:
                        print('neural gene: '+neural_gene+' tpm: '+str(neural_tpm)+'; non neural gene: '+non_neural_gene+' tpm: '+str(non_neural_tpm)+';'+' diff from neural tpm: '+str(diff))
                    non_neural_per_neural_genes_subset.append(non_neural_gene)
                    non_neural_subset_tpms[nearest_non_neural_tpm_index]=-1
            non_neural_per_neural_sites_subset_from_sites = sites_df[sites_df['bim_id'].isin(non_neural_per_neural_genes_subset)]
            print(str(len(neural_genes_subset))+' neural genes')
            print(str(len(non_neural_genes_subset))+' non_neural genes')
            print(str(len(non_neural_per_neural_genes_subset))+' non_neural_per_neural genes')
            
            print('Writing data')
            neural_outf_path = out_dir+str(x)+'fold_neural_subset'
            non_neural_outf_path = out_dir+str(x)+'fold_non_neural_subset'
            non_neural_per_neural_outf_path = out_dir+str(x)+'fold_'+str(non_neural_per_neural)+'non_neural_per_neural_subset'
            with open(neural_outf_path, 'w') as neural_fp, open(non_neural_outf_path, 'w') as non_neural_fp, open(non_neural_per_neural_outf_path, 'w') as fp:
                neural_fp.write('\n'.join(neural_genes_subset))
                non_neural_fp.write('\n'.join(non_neural_genes_subset))
                fp.write('\n'.join(non_neural_per_neural_genes_subset))    
            
            neural_sites_subset_from_sites.to_csv(neural_outf_path+'_tbl',sep='\t')
            non_neural_sites_subset_from_sites.to_csv(non_neural_outf_path+'_tbl',sep='\t')
            non_neural_per_neural_sites_subset_from_sites.to_csv(non_neural_per_neural_outf_path+'_tbl',sep='\t')
        

    if plot: # FOR PLOTING RESULTS FOR CHOCEN NEURAL/NON-NEURAL SETS


        # neural-nonneural HPM - sites types rates and editing levels distributions
        hpm_rates_dict = {}
        animals=['oct','bim','sep','squ','bob','lin']
        non_neural_hpm_res_path = 'D:/RNA_Editing_large_files_Backup_20201205/SALMON/nnpn/'
        rates_df, el_df = collect_results_for_hpm_and_plot_probs(non_neural_hpm_res_path,animals)      
        hpm_rates_dict.update({'non_neural':rates_df})    
        neural_hpm_res_path = 'D:/RNA_Editing_large_files_Backup_20201205/SALMON/neural/'
        rates_df, el_df = collect_results_for_hpm_and_plot_probs(neural_hpm_res_path,animals)
        hpm_rates_dict.update({'neural':rates_df})
    
        print('comparison of syn-div sites odds ratios (div_neural_rate/syn_neural_rate ?= div_non_neural_rate/syn_non_neural_rate):')
        for a in animals:
            a1 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'div_edited']
            a2 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'div']
            b1 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'syn_edited']
            b2 = hpm_rates_dict['neural']['strong'].loc[animals_names_dict[a],'syn']
            c1 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'div_edited']
            c2 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'div']
            d1 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'syn_edited']
            d2 = hpm_rates_dict['non_neural']['strong'].loc[animals_names_dict[a],'syn']
            z,p = calculates_odds_ratios_z_score(a1,a2,b1,b2,c1,c2,d1,d2)
            print(a+': z='+str(round(z,3))+', p='+str(round(p,3))+ ' neural '+str(round((a1/a2)/(b1/b2),4))+' nonneural '+str(round((c1/c2)/(d1/d2),4)))
        
        outfile = 'D:/RNA_Editing_large_files_Backup_20201205/SALMON/neural_and_non_neural_sites_rates_by_type.png'
        plot_odds_ratios_for_sites_groups(outfile,hpm_rates_dict, animals=animals)
     
        
# =============================================================================
#         # neural-nonneural HPM - conserved sites A>G substitution by type
#         print('comparison of C sites substitutions per type in neural/non_neural genes (chi^2):')
#         for t in ['syn','div']:
#             neural_subs = conserved_sites_subs_dict['neural'].loc[2,t+'_subs']
#             neural_sites = conserved_sites_subs_dict['neural'].loc[2,t]
#             non_neural_subs = conserved_sites_subs_dict['non_neural'].loc[2,t+'_subs']
#             non_neural_sites = conserved_sites_subs_dict['non_neural'].loc[2,t]
#             chi,p = stats.fisher_exact([[neural_subs,neural_sites],[non_neural_subs,non_neural_sites]])
#             print(t+': chi='+str(round(chi,3))+', p='+str(round(p,3)))
# =============================================================================
        

# =============================================================================
#         # neural-nonneural HPM - comparison of syn/nonsyn A>G Mutations rates ratios (for sites / background) in conserved sites for
#         general_model_dict={}
#         non_neural_path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/non_neural/4fold/coleoids_adaptive_model/'    
#         results_dfs_dict = collect_results_for_general_model_and_plot_probs(non_neural_path)
#         general_model_dict.update({'non_neural':results_dfs_dict})
#         
#         neural_path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/neural/4fold/coleoids_adaptive_model/'    
#         results_dfs_dict = collect_results_for_general_model_and_plot_probs(neural_path)
#         general_model_dict.update({'neural':results_dfs_dict})
#         
#         print('comparison sites mutations to background mutations odds ratios in neural/non_neural genes:')
#         for n in ['D','S','B']:
#             if n=='D':
#                 animals = ['sep', 'squ', 'bob', 'lin']
#             elif n=='S':
#                 animals = ['sep', 'squ']
#             elif n=='B':
#                 animals = ['bob', 'lin']
#                 
#             for a in animals:
#                 index = 'ancesNone_C_'+n+'_'+a+'_iden0.3_in_range10_liberal_average_strong_lower0.1_weak_upper0.05'
#                 a1 = general_model_dict['neural']['iden0.3_range10'].loc[index,'actual_nonsyn_eg_mutations']
#                 a2 = general_model_dict['neural']['iden0.3_range10'].loc[index,'strong_non_syn_sites']
#                 b1 = general_model_dict['neural']['iden0.3_range10'].loc[index,'non_syn_ag_mut_in_leaf']-a1
#                 b2 = general_model_dict['neural']['iden0.3_range10'].loc[index,'intermediate_non_syn_a']-a2
#                 c1 = general_model_dict['non_neural']['iden0.3_range10'].loc[index,'actual_nonsyn_eg_mutations']
#                 c2 = general_model_dict['non_neural']['iden0.3_range10'].loc[index,'strong_non_syn_sites']
#                 d1 = general_model_dict['non_neural']['iden0.3_range10'].loc[index,'non_syn_ag_mut_in_leaf']-c1
#                 d2 = general_model_dict['non_neural']['iden0.3_range10'].loc[index,'intermediate_non_syn_a']-c2
#                 z,p = calculates_odds_ratios_z_score(a1,a2,b1,b2,c1,c2,d1,d2)
#                 print('C->'+n+'->'+a+': z='+str(round(z,3))+', p='+str(round(p,3)), ' neural '+str(round((a1/a2)/(b1/b2),4))+' nonneural '+str(round((c1/c2)/(d1/d2),4)))
# =============================================================================
    

    
    
    # path_to_processed_tbl = 'D:/RNA_Editing_large_files_Backup_20201205/Phylogeny/results/Sanchez/coleoids/ST4_coleoids_orthologs_ed_sites.txt'
    # sites_expression_df = pd.read_csv(path_to_processed_tbl, sep='\t', index_col=False)
    # sites_expression_df['neural_expression'] = sites_expression_df.apply(lambda row: final_tbl_only_tpm_averaged.loc[row['O.bim_trinity_component'],'neural_average'], axis=1)
    # sites_expression_df['non_neural_expression'] = sites_expression_df.apply(lambda row: final_tbl_only_tpm_averaged.loc[row['O.bim_trinity_component'],'non_neural_average'], axis=1)
    # sites_expression_df.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/all_coleoids_sites_with_tissues_expression.xlsx')
