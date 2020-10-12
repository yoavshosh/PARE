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




def calc_average(row,tissue):

    return sum([row[c] for c in list(final_tbl_only_tpm.columns) if tissue in c])/float(len([row[c] for c in list(final_tbl_only_tpm.columns) if tissue in c]))

if __name__=='__main__':

    neural_tiss = ['optic_lobe','axial nerve_cord','subesophageal_brain','supraesophageal_brain']
    non_neural_tiss = ['ova','testes','posterior_salivary_gland','skin','suckers','embryos','retina','hepatopancreas_kidney_heart']

    
    quant_columns = ['Name','Length','EffectiveLength','TPM','NumReads']
    
    # salmon_quant_path = sys.argv[1]
    # sra_metadata_tbl = sys.argv[2]
    
    salmon_quant_path = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/SALMON/quant/'
    sra_metadata_tbl = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/SALMON/SraRunTable.txt'
    sites_tbl = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/coleoids/edited_rows_coleoids'
    above=False
    
    sra_tbl = pd.read_csv(sra_metadata_tbl,sep=',',engine='python')
    sra_tbl.set_index('Run',inplace=True)
    
    sra_names = []    
    sra_dfs = []
    quant_results=glob.glob('/'.join([salmon_quant_path,'*_quant/quant.sf']))
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
    final_tbl.to_csv('/'.join(sra_metadata_tbl.split('/')[:-1])+'/SALMON_results_all_samples.txt',sep='\t')
    final_tbl.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/SALMON_results_all_samples.xlsx')
    final_tbl_only_tpm = final_tbl[[c for c in final_tbl.columns if "TPM" in c]]
    final_tbl_only_tpm.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/SALMON_results_all_samples_tpm.xlsx')
    
    for t in sra_tbl['Tissue'].values:
        tiss = t.replace(' ','_').replace('\\','').replace(',','')
        final_tbl_only_tpm[tiss+'_average']=final_tbl_only_tpm.apply(lambda row: calc_average(row,tiss), axis=1)
    
    final_tbl_only_tpm_averaged  = final_tbl_only_tpm[[c for c in final_tbl_only_tpm if 'average' in c]]
    final_tbl_only_tpm_averaged['neural_average']=final_tbl_only_tpm_averaged.apply(lambda row: 
                                                                                    sum([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in neural_tiss])])/
                                                                                    float(len([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in neural_tiss])])), axis=1)
    
    final_tbl_only_tpm_averaged['non_neural_average']=final_tbl_only_tpm_averaged.apply(lambda row: 
                                                                                    sum([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in non_neural_tiss])])/
                                                                                    float(len([row[c] for c in list(final_tbl_only_tpm_averaged.columns) if any([t in c for t in non_neural_tiss])])), axis=1)
    
    final_tbl_only_tpm_averaged['neural_xfold'] = final_tbl_only_tpm_averaged.apply(lambda row: row['neural_average']/row['non_neural_average'] if row['non_neural_average']!=0 else np.inf, axis=1)
    final_tbl_only_tpm_averaged.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/SALMON_results_all_samples_averaged_tpm.xlsx')
    
    

    
    sites_df = pd.read_csv(sites_tbl,sep='\t', index_col=False)
    # sites_df=sites_df[np.logical_and(sites_df['nau_edited']==1,sites_df['edited_animal']==1)]
    
    
    sites_subsets = {'no_filter':len(sites_df)}
    for x in [1.5,2,2.5,3,4,5,7,10]:
        if above:
            genes_subset = final_tbl_only_tpm_averaged[np.logical_and(final_tbl_only_tpm_averaged['neural_xfold']>=x, final_tbl_only_tpm_averaged['neural_average']>1)].index
            outf_path = '/'.join(sra_metadata_tbl.split('/')[:-1])+'/'+str(x)+'fold_neural_subset'
        else:
            genes_subset = final_tbl_only_tpm_averaged[np.logical_and(final_tbl_only_tpm_averaged['neural_xfold']<x, final_tbl_only_tpm_averaged['neural_average']>1)].index
            outf_path = '/'.join(sra_metadata_tbl.split('/')[:-1])+'/'+str(x)+'fold_non_neural_subset'
        
        genes_subset = ['bim|'+i for i in genes_subset]
        sites_subset = sites_df[sites_df['bim_id'].isin(genes_subset)]
        sites_subsets.update({x:len(sites_subset)})
        with open(outf_path, 'w') as fp:
            fp.write('\n'.join(genes_subset))

    path_to_processed_tbl = 'C:/Users/shosh/OneDrive/Desktop/RNA_Editing_large_files_Backup_20200906/Phylogeny/results/Sanchez/coleoids/ST4_coleoids_orthologs_ed_sites.txt'
    sites_expression_df = pd.read_csv(path_to_processed_tbl, sep='\t', index_col=False)
    sites_expression_df['neural_expression'] = sites_expression_df.apply(lambda row: final_tbl_only_tpm_averaged.loc[row['O.bim_trinity_component'],'neural_average'], axis=1)
    sites_expression_df['non_neural_expression'] = sites_expression_df.apply(lambda row: final_tbl_only_tpm_averaged.loc[row['O.bim_trinity_component'],'non_neural_average'], axis=1)
    sites_expression_df.to_excel('/'.join(sra_metadata_tbl.split('/')[:-1])+'/all_coleoids_sites_with_tissues_expression.xlsx')
    
    