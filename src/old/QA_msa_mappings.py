# -*- coding: utf-8 -*-
"""
QA functions and operations for the ortholog mappings based on the msa results
@author: shosh
"""
import os
import glob
import pandas as pd
import numpy as np
import itertools as it
import multiprocessing as mp
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna 
import argparse
import copy


def create_unified_key(row, animals):
    
    key_str = ''
    for a in animals:
        key_str = key_str + '&&' + row[a+'_component']+';'+str(row[a+'_coding_location'])
    return key_str[2:]
    
    

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
        pair_best_hits = pair_best_hits.rename(index = str, columns = {'score':a1+'_'+a2+'_score','seq1':a1,'seq2':a2}) 
        
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


def retrive_relevant_combined_fields_for_pairs_dfs(pair,file1,file2):
    
    """
    This function reads 2 dataframes for pair - one for each full list of an animal in the pair (mapped to orthologs in the other animal)
    It returns a merged table with relevant fields per each of the animals (the returned table contain all editing sites from both animals and the mappings to each other)
    """

    relevant_fields = ['component','protein','strand','location','mm_type','aa_mm_type','editing_level','nuc','site_coding_loc_base0','site_loc_in_codon_base0','site_key']
    
    #reading dataframes and fixing combining aa before and after to one field
    df1 = pd.read_csv(file1, sep = '\t', index_col=False, low_memory=False)
    df1[pair[0]+'_aa_mm_type'] = df1.apply(lambda row: row[pair[0]+'_AA_before']+row[pair[0]+'_AA_after'],axis =1 )
    df1[pair[1]+'_aa_mm_type'] = df1.apply(lambda row: row[pair[1]+'_AA_before']+row[pair[1]+'_AA_after'],axis =1 )
    
    df2 = pd.read_csv(file2, sep = '\t', index_col=False, low_memory=False)
    df2[pair[0]+'_aa_mm_type'] = df2.apply(lambda row: row[pair[0]+'_AA_before']+row[pair[0]+'_AA_after'],axis =1 )
    df2[pair[1]+'_aa_mm_type'] = df2.apply(lambda row: row[pair[1]+'_AA_before']+row[pair[1]+'_AA_after'],axis =1 )
    
    #the relevant fields for both animals in both dataframes
    relevant_col_for_pair_df1 = [pair[0]+'_'+col for col in relevant_fields] + [pair[0]+'_'+pair[1]+'_ortholog', pair[0]+'_'+pair[1]+'_conserved', pair[0]+'_'+pair[1]+'_bitscore'] + [pair[1]+'_'+col for col in relevant_fields[1:]]
    relevant_col_for_pair_df2 = [pair[0]+'_'+col for col in relevant_fields[1:]] + [pair[1]+'_'+pair[0]+'_ortholog'] + [pair[1]+'_'+col for col in relevant_fields] 
    
    df1 = df1[relevant_col_for_pair_df1]
    df1[pair[1]+'_component'] = df1.apply(lambda row: row[pair[0]+'_'+pair[1]+'_ortholog'].split('|')[-1], axis = 1)
    df1 = df1.drop(pair[0]+'_'+pair[1]+'_ortholog', axis = 1)
    
    df2 = df2[relevant_col_for_pair_df2]
    df2[pair[0]+'_component'] = df2.apply(lambda row: row[pair[1]+'_'+pair[0]+'_ortholog'].split('|')[-1], axis = 1)
    df2 = df2.drop(pair[1]+'_'+pair[0]+'_ortholog', axis = 1)
    
    #fixing pair sites key field
    
    pair_key_col = pair[0]+'_'+pair[1]+'_coding_key'
    df1[pair[0]+'_coding_key'] = df1.apply(lambda row: str(row[pair[0]+'_component'])+';'+str(row[pair[0]+'_site_coding_loc_base0']).split('.')[0], axis=1)
    df1[pair[1]+'_coding_key'] = df1.apply(lambda row: str(row[pair[1]+'_component'])+';'+str(row[pair[1]+'_site_coding_loc_base0']).split('.')[0], axis=1)
    df2[pair[0]+'_coding_key'] = df2.apply(lambda row: str(row[pair[0]+'_component'])+';'+str(row[pair[0]+'_site_coding_loc_base0']).split('.')[0], axis=1)
    df2[pair[1]+'_coding_key'] = df2.apply(lambda row: str(row[pair[1]+'_component'])+';'+str(row[pair[1]+'_site_coding_loc_base0']).split('.')[0], axis=1)
    df1[pair_key_col] = df1.apply(lambda row: str(row[pair[0]+'_coding_key'])+'|'+str(row[pair[1]+'_coding_key']), axis=1)
    df2[pair_key_col] = df2.apply(lambda row: str(row[pair[0]+'_coding_key'])+'|'+str(row[pair[1]+'_coding_key']), axis=1)
    
    return pd.concat([df1,df2[~df2[pair_key_col].isin(list(df1[pair_key_col]))]], sort=False, ignore_index = True)


def compair_tables(pairwise_mappings, msa_mappings, super_ortholog_genes, pair):
    
    """
    compare the pairwise alignment results to the msa results for rach pair of animals
    return the sites which are not within the overlap
    """
        
    def create_new_pairwise_pair_key_for_comparison(row,a1,a2):
        if pd.isnull(row[a1+'_site_coding_loc_base0']):
            a1_key = row[a1+'_component']+';'+'gap'
        else:
            a1_key = row[a1+'_component']+';'+str(round(row[a1+'_site_coding_loc_base0'],0))
        a1_key = a1_key.split('.')[0]
        
        if pd.isnull(row[a2+'_site_coding_loc_base0']):
            a2_key = row[a2+'_component']+';'+'gap'
        else:
            a2_key = row[a2+'_component']+';'+str(round(row[a2+'_site_coding_loc_base0'],0))
        a2_key = a2_key.split('.')[0]
    
        return a1_key+'|'+a2_key
        
    a1 = pair[0]
    a2 = pair[1]
 
    pairwise_mappings[a1+'_component'] = pairwise_mappings.apply(lambda row: a1+'|'+row[a1+'_component'],axis=1)
    pairwise_mappings[a2+'_component'] = pairwise_mappings.apply(lambda row: a2+'|'+row[a2+'_component'],axis=1)
    pairwise_mappings = pairwise_mappings[np.logical_or(pairwise_mappings[a1+'_component'].isin(list(super_ortholog_genes[a1])),pairwise_mappings[a2+'_component'].isin(list(super_ortholog_genes[a2])))]
    pairwise_mappings[a1+'_'+a2+'_key'] = pairwise_mappings.apply(lambda row: create_new_pairwise_pair_key_for_comparison(row,a1,a2), axis=1)

#    msa_mappings = msa_mappings[np.logical_or(msa_mappings[a1+'_coding_location']!='gap',msa_mappings[a2+'_coding_location']!='gap')]
    msa_mappings = msa_mappings[np.logical_or(msa_mappings[a1+'_edited'],msa_mappings[a2+'_edited'])]
    msa_mappings[a1+'_'+a2+'_key'] = msa_mappings.apply(lambda row: row[a1+'_component']+';'+str(row[a1+'_coding_location'])+'|'+row[a2+'_component']+';'+str(row[a2+'_coding_location']), axis=1)
    
    msa_mapping_not_overlap = msa_mappings[~msa_mappings[a1+'_'+a2+'_key'].isin(list(pairwise_mappings[a1+'_'+a2+'_key']))]
    pairwise_mappings_not_overlap = pairwise_mappings[~pairwise_mappings[a1+'_'+a2+'_key'].isin(list(msa_mappings[a1+'_'+a2+'_key']))]
    
    print(str(pair))
    print('pairwise super ortholog mappings (where at least one animal is edited): ' + str(len(pairwise_mappings)))
    print('pairwise not in msa: ' + str(len(pairwise_mappings_not_overlap)))
    print('msa super ortholog mappings (where at least one animal is edited): ' + str(len(msa_mappings)))
    print('msa not in pairwise: ' + str(len(msa_mapping_not_overlap)))
    print('\n')
    
    return {pair:{'msa_mapping_not_overlap':msa_mapping_not_overlap,'pairwise_mappings_not_overlap':pairwise_mappings_not_overlap}}
    


if __name__ == '__main__':
    
    pairwise_path = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/results_from_blast_parsing_using_my_scripts/'
#    grand_locations_table_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/results/all_sites_mappings.txt'
    grand_locations_proteins_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/all_sites_mappings_mrna_new.txt'
    grand_locations_mrna_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/all_sites_mappings_proteins_new.txt'
    all_best_hits_file = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/all_best_hits.txt'
    
    animals = ['oct','squ','bim','nau','sep','bob','lin']

    print('Reading 2-way best hits file for animals: ' + str(animals))
    all_best_hits = pd.read_csv(all_best_hits_file, sep = '\t', names = ['seq1','seq2','score'])
    all_best_hits['animal_1'] = all_best_hits.apply(lambda row: row['seq1'][:3], axis=1)
    all_best_hits['animal_2'] = all_best_hits.apply(lambda row: row['seq2'][:3], axis=1)
    super_ortholog_genes = create_super_orthologs_df(all_best_hits,animals)
    print('Super ortholog genes: ' + str(len(super_ortholog_genes)))

    print('Collecting pairs pairwise mapping data')    
    pairwise_pairs_dict = {}
    for pair in it.combinations(animals,2):
        file1 = pairwise_path + pair[0]+'_editing_sites_'+pair[1]+'_orthologs.txt'
        file2 = pairwise_path + pair[1]+'_editing_sites_'+pair[0]+'_orthologs.txt'
        pair_complete_df = retrive_relevant_combined_fields_for_pairs_dfs(pair,file1,file2)
        pairwise_pairs_dict.update({pair[0]+'_'+pair[1]:pair_complete_df})
        
    print('Reading MSA mappings table for all animals')
    grand_locations_table_proteins = pd.read_csv(grand_locations_proteins_file, sep = '\t')
    grand_locations_table_mrna = pd.read_csv(grand_locations_mrna_file, sep = '\t')
    
    
    pairwise_vs_mrna_msa_notoverlap_dict = {}
    pairwise_vs_proteins_msa_notoverlap_dict = {}
    pairwise_pairs_dict_2 = copy.deepcopy(pairwise_pairs_dict)
    pairwise_pairs_dict_3 = copy.deepcopy(pairwise_pairs_dict)
    for pair in it.combinations(animals,2):
        
        try:
            pairwise_mappings_2 = pairwise_pairs_dict_2[pair[0]+'_'+pair[1]]
            pairwise_mappings_3 = pairwise_pairs_dict_3[pair[0]+'_'+pair[1]]
        except KeyError:
            pairwise_mappings_2 = pairwise_pairs_dict_2[pair[1]+'_'+pair[0]]
            pairwise_mappings_3 = pairwise_pairs_dict_3[pair[1]+'_'+pair[0]]
            
        print('pairwise vs mrna msa')
        pair_notoverlap_dict = compair_tables(pairwise_mappings_2,grand_locations_table_mrna,super_ortholog_genes,pair)
        pairwise_vs_mrna_msa_notoverlap_dict.update(pair_notoverlap_dict)
        print('pairwise vs proteins msa')
        pair_notoverlap_dict = compair_tables(pairwise_mappings_3,grand_locations_table_proteins,super_ortholog_genes,pair)
        pairwise_vs_proteins_msa_notoverlap_dict.update(pair_notoverlap_dict)
        
        
        
grand_locations_table_mrna['unified_key'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, animals), axis=1)        
grand_locations_table_proteins['unified_key'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, animals), axis=1)
mrna = list(grand_locations_table_mrna['unified_key'])
prot = list(grand_locations_table_proteins['unified_key'])
mrna_not_in_prot = [k for k in mrna if k not in prot]
prot_not_in_mrna = [k for k in prot if k not in mrna]


grand_locations_table_mrna['unified_key_wo_nau'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','bob','lin']), axis=1)        
grand_locations_table_proteins['unified_key_wo_nau'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','bob','lin']), axis=1)
mrna_wo_nau = list(grand_locations_table_mrna['unified_key_wo_nau'])
prot_wo_nau = list(grand_locations_table_proteins['unified_key_wo_nau'])
mrna_not_in_prot_wo_nau = [k for k in mrna_wo_nau if k not in prot_wo_nau]
prot_not_in_mrna_wo_nau = [k for k in prot_wo_nau if k not in mrna_wo_nau]


grand_locations_table_mrna['unified_key_wo_bob'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','nau','lin']), axis=1)        
grand_locations_table_proteins['unified_key_wo_bob'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','nau','lin']), axis=1)
mrna_wo_bob = list(grand_locations_table_mrna['unified_key_wo_bob'])
prot_wo_bob = list(grand_locations_table_proteins['unified_key_wo_bob'])
mrna_not_in_prot_wo_bob = [k for k in mrna_wo_bob if k not in prot_wo_bob]
prot_not_in_mrna_wo_bob = [k for k in prot_wo_bob if k not in mrna_wo_bob]


grand_locations_table_mrna['unified_key_wo_lin'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','nau','bob']), axis=1)        
grand_locations_table_proteins['unified_key_wo_lin'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['oct','squ','bim','sep','nau','bob']), axis=1)
mrna_wo_lin = list(grand_locations_table_mrna['unified_key_wo_lin'])
prot_wo_lin = list(grand_locations_table_proteins['unified_key_wo_lin'])
mrna_not_in_prot_wo_lin = [k for k in mrna_wo_lin if k not in prot_wo_lin]
prot_not_in_mrna_wo_lin = [k for k in prot_wo_lin if k not in mrna_wo_lin]


grand_locations_table_mrna['unified_key_wo_sep'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['oct','squ','bim','lin','nau','bob']), axis=1)        
grand_locations_table_proteins['unified_key_wo_sep'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['oct','squ','bim','lin','nau','bob']), axis=1)
mrna_wo_sep = list(grand_locations_table_mrna['unified_key_wo_sep'])
prot_wo_sep = list(grand_locations_table_proteins['unified_key_wo_sep'])
mrna_not_in_prot_wo_sep = [k for k in mrna_wo_sep if k not in prot_wo_sep]
prot_not_in_mrna_wo_sep = [k for k in prot_wo_sep if k not in mrna_wo_sep]


grand_locations_table_mrna['unified_key_wo_bim'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['oct','squ','sep','lin','nau','bob']), axis=1)        
grand_locations_table_proteins['unified_key_wo_bim'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['oct','squ','sep','lin','nau','bob']), axis=1)
mrna_wo_bim = list(grand_locations_table_mrna['unified_key_wo_bim'])
prot_wo_bim = list(grand_locations_table_proteins['unified_key_wo_bim'])
mrna_not_in_prot_wo_bim = [k for k in mrna_wo_bim if k not in prot_wo_bim]
prot_not_in_mrna_wo_bim = [k for k in prot_wo_bim if k not in mrna_wo_bim]


grand_locations_table_mrna['unified_key_wo_oct'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['bim','squ','sep','lin','nau','bob']), axis=1)        
grand_locations_table_proteins['unified_key_wo_oct'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['bim','squ','sep','lin','nau','bob']), axis=1)
mrna_wo_oct = list(grand_locations_table_mrna['unified_key_wo_oct'])
prot_wo_oct = list(grand_locations_table_proteins['unified_key_wo_oct'])
mrna_not_in_prot_wo_oct = [k for k in mrna_wo_oct if k not in prot_wo_oct]
prot_not_in_mrna_wo_oct = [k for k in prot_wo_oct if k not in mrna_wo_oct]


grand_locations_table_mrna['unified_key_wo_squ'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['bim','oct','sep','lin','nau','bob']), axis=1)        
grand_locations_table_proteins['unified_key_wo_squ'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['bim','oct','sep','lin','nau','bob']), axis=1)
mrna_wo_squ = list(grand_locations_table_mrna['unified_key_wo_squ'])
prot_wo_squ = list(grand_locations_table_proteins['unified_key_wo_squ'])
mrna_not_in_prot_wo_squ = [k for k in mrna_wo_squ if k not in prot_wo_squ]
prot_not_in_mrna_wo_squ = [k for k in prot_wo_squ if k not in mrna_wo_squ]


grand_locations_table_mrna['unified_key_wo_bob_lin'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['bim','oct','sep','nau','squ']), axis=1)        
grand_locations_table_proteins['unified_key_wo_bob_lin'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['bim','oct','sep','nau','squ']), axis=1)
mrna_wo_bob_lin = list(grand_locations_table_mrna['unified_key_wo_bob_lin'])
prot_wo_bob_lin = list(grand_locations_table_proteins['unified_key_wo_bob_lin'])
mrna_not_in_prot_wo_bob_lin = [k for k in mrna_wo_bob_lin if k not in prot_wo_bob_lin]
prot_not_in_mrna_wo_bob_lin = [k for k in prot_wo_bob_lin if k not in mrna_wo_bob_lin]


grand_locations_table_mrna['unified_key_wo_bob_lin_squ'] = grand_locations_table_mrna.apply(lambda row: create_unified_key(row, ['bim','oct','sep','nau']), axis=1)        
grand_locations_table_proteins['unified_key_wo_bob_lin_squ'] = grand_locations_table_proteins.apply(lambda row: create_unified_key(row, ['bim','oct','sep','nau']), axis=1)
mrna_wo_bob_lin_squ = list(grand_locations_table_mrna['unified_key_wo_bob_lin_squ'])
prot_wo_bob_lin_squ = list(grand_locations_table_proteins['unified_key_wo_bob_lin_squ'])
mrna_not_in_prot_wo_bob_lin_squ = [k for k in mrna_wo_bob_lin_squ if k not in prot_wo_bob_lin_squ]
prot_not_in_mrna_wo_bob_lin_squ = [k for k in prot_wo_bob_lin_squ if k not in mrna_wo_bob_lin_squ]


"""
Tree based on file msa_results_for_super_orthologs_3199 (nucleotides)
((OCT|COMP168018_C0_SEQ1:0.002074,BIM|COMP99991_C0_SEQ1:0.005434)N2:0.335050,(SEP|COMP257715_C0_SEQ1:0.064550,(SQU|COMP140236_C0_SEQ1:0.055089,(BOB|COMP374597_C2_SEQ1:0.120148,LIN|COMP766105_C0_SEQ2:0.076994)N5:0.026150)N4:0.022458)N3:0.181966,NAU|COMP163206_C0_SEQ1:0.571687)N1;
"""