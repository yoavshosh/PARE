#import re
#import os
import sys
import pandas as pd
import numpy as np
import itertools as it
import multiprocessing as mp
#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#import statsmodels.stats.multitest as p_adjust
#from pylab import text
#from scipy import stats
#from collections import deque
#from functools import reduce
#from matplotlib import colors
#from matplotlib.colors import LogNorm
#from heapq import nsmallest
#from Bio import SeqIO
#import xlsxwriter
import pickle

def create_ambigouos_mappings_dict_parallel(animals,pairs_dict):
    
    output = mp.Queue()
    processes = []
    ambigouos_mappings = {}
    trios_positions = {}
    pos=0
    for pair in pairs_dict.keys():
        p = pair.split('_')
        pair_df = pairs_dict[p[0]+'_'+p[1]]
        for a in [i for i in animals if i not in p]:
            pos+=1
            try:
                first_with_third = pairs_dict[p[0]+'_'+a]
            except KeyError:
                first_with_third = pairs_dict[a+'_'+p[0]]
            try:
                second_with_third = pairs_dict[p[1]+'_'+a]
            except KeyError:
                second_with_third = pairs_dict[a+'_'+p[1]]
            
            animals = p+[a]
            processes.append(mp.Process(target=check_ambigouos_mapping_in_third_animal, args=(animals, pair_df, first_with_third, second_with_third, pos, output)))
            trios_positions.update({pos:'_'.join(animals)})
            
    for p in processes:
        p.start()
    for p in processes:
        p.join()
        
    results = [output.get() for p in processes]
    
    for res in results:
        ambigouos_mappings.update({trios_positions[res[0]]:res[1]})
        
    return ambigouos_mappings
    

def check_ambigouos_mapping_in_third_animal(animals, pair_df, first_with_third, second_with_third, pos, output):
    """
    animal = a list of 3 animals
    check for each orthologs in pair_df of first 2 animals if sites are mapped to different locations in third animal
    This is rather for QA and not used to build grand unified table
    """
    a1 = animals[0]
    a2 = animals[1]
    a3 = animals[2]
    
    print(a1+'_'+a2 +': checking mappings to ' + a3)
    
    connected_rows = pair_df[np.logical_and(pair_df[a1+'_component']!='no_ortholog',pair_df[a2+'_component']!='no_ortholog')]
    ambigouos_mapping_to_third_animal = {}
    for index, row in connected_rows.iterrows():
        try:
#            if row[a1+'_component'] != 'no_ortholog' and row[a2+'_component'] != 'no_ortholog':
           if 'nan' not in row[a1+'_coding_key'] + row[a2+'_coding_key']:
               a1_row_in_third = first_with_third[first_with_third[a1+'_coding_key']==row[a1+'_coding_key']].squeeze()
               a2_row_in_third = second_with_third[second_with_third[a2+'_coding_key']==row[a2+'_coding_key']].squeeze()
              
               if len(a1_row_in_third) and len(a2_row_in_third):
    
                   a1_a3_coding_ortholog = a1_row_in_third[a3+'_coding_key']
                   a2_a3_coding_ortholog = a2_row_in_third[a3+'_coding_key']
                   if a1_a3_coding_ortholog != a2_a3_coding_ortholog:    
                       ambigouos_mapping_to_third_animal.update({a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key']:(a1_a3_coding_ortholog,a2_a3_coding_ortholog)})
#                       print(a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key'] +' does not have a uniqe connection to ' + a3)
                   else:
#                       print(a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key'] +' have a uniqe connection to ' + a3+';'+a2_a3_coding_ortholog)                      
                       pass
           else:
               a1_a3_ortholog = first_with_third[first_with_third[a1+'_component']==row[a1+'_component']]
               a2_a3_ortholog = second_with_third[second_with_third[a2+'_component']==row[a2+'_component']]
               
               if len(a1_a3_ortholog) and len(a2_a3_ortholog):
            
                   a1_a3_ortholog_comp = a1_a3_ortholog.iloc[0,:][a3+'_component']
                   a2_a3_ortholog_comp = a2_a3_ortholog.iloc[0,:][a3+'_component']
                  
                   if a1_a3_ortholog_comp != a2_a3_ortholog_comp:
                       ambigouos_mapping_to_third_animal.update({a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key']:(a1_a3_ortholog_comp,a2_a3_ortholog_comp)})
#                       print(a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key'] +' does not have a uniqe connection to ' + a3)
                   else:
#                       print(a1+';'+row[a1+'_coding_key']+'|'+a2+';'+row[a2+'_coding_key'] +' have a uniqe connection to ' + a3+';'+a2_a3_ortholog_comp)
                       pass
        except:
            print('Error on row ' + str(index))
            break
              
    output.put((pos,ambigouos_mapping_to_third_animal))
#    return ambigouos_mapping_to_third_animal
    

def retrive_relevant_combined_fields_for_pairs_dfs(pair,file1,file2,reading_output):
    
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

#    pair_df = pd.concat([df1,df2[~df2[pair_key_col].isin(list(df1[pair_key_col]))]], sort=False, ignore_index = True)
#    reading_output.put((pair,pair_df))

    
def create_full_table_from_dict(pairs_dict,animals):
    
    relevant_fields = ['component','protein','strand','location','mm_type','aa_mm_type','editing_level','nuc','site_coding_loc_base0','site_loc_in_codon_base0','site_key']
    
    conserved_dict = {}
    df_all_animals = None
    for pair, df in pairs_dict.items():
        
        a1 = pair.split('_')[0]
        a2 = pair.split('_')[1]
        
        if all(a in animals for a in [a1,a2]):
        
            df = df[df[a1+'_'+a2+'_conserved']==True]       
            df[a1+'_'+a2+'_recoding_event'] = df.apply(lambda row: row[a1+'_aa_mm_type'][0]!=row[a1+'_aa_mm_type'][1] and row[a2+'_aa_mm_type'][0]!=row[a2+'_aa_mm_type'][1], axis = 1)
            conserved_dict.update({pair:df})
            
            if df_all_animals is None:
                df_all_animals = df
            else:
                #both animals already exist in df_all_animals - merge by both animals kyes and eventually take from temp df the alignment_identity column only
                if all(a+'_coding_key' in df_all_animals.columns for a in [a1,a2]):
                    print('both animals are already in df')
                    columns_to_use = [a1 + '_coding_key',a2 + '_coding_key'] + [a1+'_'+a2+'_bitscore',a1+'_'+a2+'_recoding_event']
                    df_all_animals = df_all_animals.merge(df[columns_to_use], on = [a1 + '_coding_key',a2 + '_coding_key'], how = 'inner')
                #only animal1 exits in df_all_animals - merge by animal1 key and eventually take from temp df all animal2 columns (and alignment identity)
                elif a1 + '_coding_key' in df_all_animals.columns and a2 + '_coding_key' not in df_all_animals.columns:
                    print('animal_1 (' + a1 + ') is already in df')
                    columns_to_use = [a2 + '_' + x for x in relevant_fields] + [a1 + '_coding_key',a2 + '_coding_key'] + [a1+'_'+a2+'_bitscore',a1+'_'+a2+'_recoding_event']
                    df_all_animals = df_all_animals.merge(df[columns_to_use], on = a1 + '_coding_key', how = 'inner')
                #same logic as last condition but now only animal2 is already exist in df_all_animals
                elif a2 + '_coding_key' in df_all_animals.columns and a1 + '_coding_key' not in df_all_animals.columns:
                    print('animal_2 (' + a2 + ') is already in df')
                    columns_to_use = [a1 + '_' + x for x in relevant_fields] + [a1 + '_coding_key',a2 + '_coding_key'] + [a1+'_'+a2+'_bitscore',a1+'_'+a2+'_recoding_event']
                    df_all_animals = df_all_animals.merge(df[columns_to_use], on = a2 + '_coding_key', how = 'inner')
                if all(a+'_coding_key' not in df_all_animals.columns for a in [a1,a2]):
                    print('NO INTERSECTION OF ANIMALS')
                    
    conserved_dict.update({'conserved_across_all':df_all_animals})
    return conserved_dict



# =============================================================================
# def fined_pair_mappings_in_other_animals(pairs_dict):
#     """
#     This function iterate over dfs in pairs_dict (created by retrive_relevant_combined_fields_for_pairs_dfs)
#     for each pair_df it mapps each of the two fields containing coding kyes for thetwo animals, to all other orthologs in other animals using the other dfs from pairs_dict
#     """
#     
#     def find_mapping_to_animal(row, a1, other_animal, a1_with_other_animal):
#         
#         possible_mappings = a1_with_other_animal[a1_with_other_animal[a1+'_coding_key']==row[a1+'_coding_key']]
#         possible_coding_keys = possible_mappings[[other_animal+'_coding_key']]
#         if len(set(list(possible_coding_keys)))>1:
#             raise Exception(a1 + '|' + row[a1+'_coding_key'] + ' has an unambigouos mapping to ' + other_animal)
#         else:
#             row[a1+'_'+other_animal+'_ortholog_coding_key'] = list(possible_coding_keys)[0]
#         
#         return row
#     
#     def find_animals_mappings_in_other_pairs(a1, a2, a1_df, pairs_dict):
#         
#         for pair, df in pairs_dict.items():
#             animals = pair.split('_')
#             if a1 in animals and a2 not in animals:
#                 other_animal = animals[animals.index(a)-1]
#                 print('mapping ' + a1 + ' to ' + other_animal)
#                 a1_df = a1_df.apply(lambda row: find_mapping_to_animal(row, a1, other_animal, df), axis=1)
#         
#         return a1_df
#                 
#     new_dict = {}
#     for pair, pairs_df in pairs_dict.items():
#         print('mapping ' + pair + ' to other animals')
#         animals = pair.split('_')
#         new_df = pairs_df.copy()
#         for a in animals:
#             new_df = find_animals_mappings_in_other_pairs(a, animals[animals.index(a)-1], new_df, pairs_dict)
#         new_dict.update({pair:new_df})
#         
#     return new_dict
# =============================================================================
    
        
def fined_pair_mappings_in_other_animals_parallel(pairs_dict):
    """
    This function iterate over dfs in pairs_dict (created by retrive_relevant_combined_fields_for_pairs_dfs)
    for each pair_df it mapps each of the two fields containing coding kyes for thetwo animals, to all other orthologs in other animals using the other dfs from pairs_dict
    """
    
    def find_mapping_to_animal(row, a1, other_animal, a1_with_other_animal):
        """
        for row with a1 animal, return mapping to other_animal from a1_with_other_animal
        """
        possible_mappings = a1_with_other_animal[a1_with_other_animal[a1+'_coding_key']==row[a1+'_coding_key']]
        possible_coding_keys = possible_mappings[[other_animal+'_coding_key']]
        if len(set(list(possible_coding_keys)))>1:
            raise Exception(a1 + '|' + row[a1+'_coding_key'] + ' has an unambigouos mapping to ' + other_animal)
        else:
            row[a1+'_'+other_animal+'_mapping'] = list(possible_coding_keys)[0]
        return row
    
    def find_animals_mappings_in_other_pairs(a1, a2, a1_df, pairs_dict, pos, output):
        """
        return a1_a2 pair dataframe with all mappings of a1 to other animals
        """
        for pair, df in pairs_dict.items():
            animals = pair.split('_')
            if a1 in animals and a2 not in animals:
                other_animal = animals[animals.index(a1)-1]
                print(a1+'_'+a2 + ': mapping ' + a1 + ' to ' + other_animal)
                a1_df = a1_df.apply(lambda row: find_mapping_to_animal(row, a1, other_animal, df), axis=1)
        output.put((pos,a1_df))
       
    #for each pair crate 2 processes. one for each animal each process mapps the animal to all other animal and returns the mapped dataframe
    output = mp.Queue()
    processes = []
    for pos, item in enumerate(pairs_dict.items()):
        pair = item[0]
        df = item[1]
        animals = pair.split('_')
        for a in animals:
            processes.append(mp.Process(target=find_animals_mappings_in_other_pairs, args=(a, animals[animals.index(a)-1], df, pairs_dict, pos, output)))
    
    
    for p in processes:
        p.start()
        
    for p in processes:
        p.join()
    
    results = [output.get() for p in processes]
    
    return results
        
def create_new_pairs_dict_from_other_animal_mappings_results(mapped_dfs, pairs_keys):
    n_pairs = len(mapped_dfs)/2
    
    new_dict = {}
    for i, key in pairs_keys:
        dfs_to_connect = []
        for pos, df in mapped_dfs:
            if pos == i:
                dfs_to_connect.append(df)
        if len(dfs_to_connect)!=2:
            raise Exception(key + ' dose not have 2 dataframes')
        else:
            columns = [col for col in dfs_to_connect[1].columns if 'mapping' in col] + [key + '_coding_key']
            cosolidated_df = dfs_to_connect[0].merge(dfs_to_connect[1][columns], on = key+'_coding_key', how = 'outer')
            new_dict.update({key:cosolidated_df})
            if any(len(cosolidated_df)!=len(dfs_to_connect[0]),len(cosolidated_df)!=len(dfs_to_connect[1]),len(dfs_to_connect[0])!=len(dfs_to_connect[1])):
                raise Exception(key + ' lengths of dataframes are not equal')
            
    return new_dict
        
    

    
def connect_all_dfs_containing_specific_animal(animal,pairs_dict):
    
    connected_df = None
    for pair, df in pairs_dict.items():
        if animal in pair:
            if connected_df is None:
                connected_df = df
            else:
                pair_list = pair.split('_')
                if pair_list[0] == animal:
                    joined_animal = pair_list[1]
                else:
                    joined_animal = pair_list[0]
                
                df_only_exact_locations_in_animal = df[np.logical_and(df[animal+'_component']!='no_ortholog',df[animal+'_site_coding_loc_base0'].notnull())]
                
                columns_to_merge = ['component','protein','strand','location','mm_type','aa_mm_type','editing_level','nuc','site_coding_loc_base0','site_loc_in_codon_base0','site_key','coding_key']
                columns_to_merge = [joined_animal+'_'+col for col in columns_to_merge] + [animal+'_coding_key']
                connected_df = connected_df.merge(df_only_exact_locations_in_animal[columns_to_merge], on = animal+'_coding_key', how = 'outer')
    
    return connected_df
    
def create_mappings_dict(animals, pairs_dict):
    
    grand_dict = {}
    for a in animals:
        print('checking mappings for ' + a)
        pairs = [key for key in pairs_dict.keys() if a in key]
        for p in pairs:
            animals_pair = p.split('_')
            other_animal = animals_pair[animals_pair.index(a)-1]
            print(other_animal)
            all_sites_for_animal = pairs_dict[p][[a+'_site_key',a+'_coding_key',other_animal+'_coding_key']]  
            all_sites_for_animal = all_sites_for_animal[~all_sites_for_animal[a+'_site_key'].isnull()]
            grand_dict.update({a+'_'+other_animal:dict(zip(all_sites_for_animal[a+'_coding_key'], all_sites_for_animal[other_animal+'_coding_key']))})        

    return grand_dict


if __name__ == '__main__':
    
#    files_path = sys.argv[1]
#    files_path = 'E:/RNA_editing_Large_files/orthomcl/compliantFasta_noa_with_lin/results/'
    files_path = 'E:/RNA_editing_Large_files/orthomcl/orthomcl_7_species/results_from_blast_parsing_using_my_scripts/'
    
    animals = ['oct','sep','squ','bim','lin','bob','nau']
    conserved_animals = ['oct','bim','sep','squ','bob','lin']

    print('Collecting pairs data')    
    reading_output = mp.Queue()
    reading_processes = []
    pairwise_pairs_dict = {}
    for pair in it.combinations(animals,2):
        file1 = files_path + pair[0]+'_editing_sites_'+pair[1]+'_orthologs.txt'
        file2 = files_path + pair[1]+'_editing_sites_'+pair[0]+'_orthologs.txt'
        pair_complete_df = retrive_relevant_combined_fields_for_pairs_dfs(pair,file1,file2,reading_output)
        pairwise_pairs_dict.update({pair[0]+'_'+pair[1]:pair_complete_df})
#        reading_processes.append(mp.Process(target=retrive_relevant_combined_fields_for_pairs_dfs, args=(pair, file1, file2, reading_output)))


#    for p in reading_processes:
#        p.start()
#    for p in reading_processes:
#        p.join()
#    reading_results = [reading_output.get() for p in reading_processes]
#    for res in reading_results:
#        pairwiswe_pairs_dict.update({res[0]:res[1]})
#        
#    print('Saving pairs dictionary')
#    np.save(files_path+'pairwiswe_pairs_dict.npy', pairwiswe_pairs_dict)
#    
#    conserved_dict = create_full_table_from_dict(pairwiswe_pairs_dict,conserved_animals)
#    print('Saving conserved dictionary')
#    np.save(files_path+'conserved_dict.npy', conserved_dict)
#    excel_writer = pd.ExcelWriter(files_path+'_'.join(conserved_animals)+'.xlsx',engine='xlsxwriter')
#    for dataframe, sheet in zip(conserved_dict.values(),conserved_dict.keys()):
#        dataframe.to_excel(excel_writer, sheet_name=sheet, startrow=0 , startcol=0, index=False) 
#    
#    pairwiswe_pairs_dict = np.load(files_path+'pairwiswe_pairs_dict.npy', encoding = 'latin1')  

#    print('Checking ambigouos mapping of pairs to other animals')
#    ambigouos_mappings = create_ambigouos_mappings_dict_parallel(animals,pairwiswe_pairs_dict)
#    print('Saving ambigouos mapping to pickle')
#    np.save(files_path+'ambigouos_mappings.npy', ambigouos_mappings)
    
    
#    print('Finding all mappings for each animal in each pair') 
#    dataframes_with_other_animals_mappings = fined_pair_mappings_in_other_animals_parallel(pairwiswe_pairs_dict)
#    print('Merging pairs results to new pairs dictionary')
#    new_pairs_dict = create_new_pairs_dict_from_other_animal_mappings_results(dataframes_with_other_animals_mappings, list(pairwiswe_pairs_dict.keys()))
#    np.save('new_pairs_dict.npy', new_pairs_dict)
#
#    with open(files_path+ 'pairs_results_all_mappings.pkl', 'wb') as f:
#        pickle.dump(results, f)
#    np.save(files_path+'pairs_results_all_mappings.npy', results)
        
    

#    print('Saving dataframes dictionary')  
#    with open(files_path+ 'pairwiswe_pairs_dict.pkl', 'wb') as f:
#        pickle.dump(pairwiswe_pairs_dict, f)
        

#        
#    with open(files_path+ 'pairwiswe_pairs_dict.pkl', 'rb') as handle:
#        pairwiswe_pairs_dict = pickle.load(handle)
    
#
#    ambigouos_mappings = {}
#    for pair in pairwiswe_pairs_dict.keys():
#        p = pair.split('_')
#        pair_df = pairwiswe_pairs_dict[p[0]+'_'+p[1]]
#        for a in [i for i in animals if i not in p]:
#            try:
#                first_with_third = pairwiswe_pairs_dict[p[0]+'_'+a]
#            except KeyError:
#                first_with_third = pairwiswe_pairs_dict[a+'_'+p[0]]
#            try:
#                second_with_third = pairwiswe_pairs_dict[p[1]+'_'+a]
#            except KeyError:
#                second_with_third = pairwiswe_pairs_dict[a+'_'+p[1]]
#            
#            if pair+'_'+a not in ambigouos_mappings.keys():    
#                print(pair+'_'+a)
#                ambigouos_mapping_in_trio = check_ambigouos_mapping_in_third_animal(p+[a],pair_df,first_with_third,second_with_third)
#                ambigouos_mappings.update({pair+'_'+a:ambigouos_mapping_in_trio})
#            
#            
#    with open(files_path+ 'ambigouos_mappings.pkl', 'wb') as f:
#        pickle.dump(ambigouos_mappings, f)
#        
#    with open(files_path+ 'ambigouos_mappings.pkl', 'rb') as handle:
#        ambigouos_mappings = pickle.load(handle)
                
#    ambigouos_mappings_double_mapped_only = {}
#    for trio, trio_dict in ambigouos_mappings.items():
#        sub_dict = {}
#        for k,v in trio_dict.items():
#            if 'no_ortholog;nan' not in v and 'no_ortholog' not in v :
#                sub_dict.update({k:v})
#        ambigouos_mappings_double_mapped_only.update({trio:sub_dict})
#                
#    ambigouos_mappings_double_mapped_exact_locations_only = {} 
#    for trio, trio_dict in ambigouos_mappings_double_mapped_only.items():
#        sub_dict = {}
#        for k,v in trio_dict.items():
#            if 'nan' not in k:
#                sub_dict.update({k:v})
#        ambigouos_mappings_double_mapped_exact_locations_only.update({trio:sub_dict})
#        
#    for trio, trio_dict in ambigouos_mappings_double_mapped_exact_locations_only.items():
#        n_same = 0
#        n_diff = 0
#        for k,v in trio_dict.items():
#            if v[0].split(';')[0] == v[1].split(';')[0]:
#                n_same+=1
#            else:
#                n_diff+=1
#        print(trio + ' n_same:'+str(n_same) + 'n_diff:'+str(n_diff))
#        
#    oct_sep = pairwiswe_pairs_dict['oct_sep']
#    
#    oct_sep = oct_sep[['nan' not in x for x in list(oct_sep['oct_coding_key'])]]
#    oct_sep = oct_sep[['nan' not in x for x in list(oct_sep['sep_coding_key'])]] 
#    len(set(oct_sep['oct_sep_coding_key']))
    
    

                
    
        
        
        
        
        
        
        
#columns_from_input_dfs = ['squ_component',
# 'squ_protein',
# 'squ_location',
# 'squ_mm_type',
# 'squ_DNA_A',
# 'squ_DNA_T',
# 'squ_DNA_G',
# 'squ_DNA_C',
# 'squ_RNA_A',
# 'squ_RNA_T',
# 'squ_RNA_G',
# 'squ_RNA_C',
# 'squ_Trinity',
# 'squ_RNA_coverage',
# 'squ_DNA_coverage',
# 'squ_p_val',
# 'squ_AA_before',
# 'squ_AA_after',
# 'squ_type',
# 'squ_protein_length',
# 'squ_editing_level',
# 'squ_strand',
# 'squ_trinity_exists',
# 'squ_nuc',
# 'squ_location_in_area',
# 'squ_site_area',
# 'squ_site_coding_loc_base0',
# 'squ_site_loc_in_codon_base0',
# 'squ_original_codon',
# 'squ_site_key',
# 'squ_bim_ortholog',
# 'squ_bim_alignment_identity',
# 'squ_bim_ortholog_coding_location',
# 'squ_bim_ortholog_site_key',
# 'bim_component',
# 'bim_protein',
# 'bim_location',
# 'bim_mm_type',
# 'bim_DNA_A',
# 'bim_DNA_T',
# 'bim_DNA_G',
# 'bim_DNA_C',
# 'bim_RNA_A',
# 'bim_RNA_T',
# 'bim_RNA_G',
# 'bim_RNA_C',
# 'bim_Trinity',
# 'bim_RNA_coverage',
# 'bim_DNA_coverage',
# 'bim_p_val',
# 'bim_AA_before',
# 'bim_AA_after',
# 'bim_type',
# 'bim_protein_length',
# 'bim_editing_level',
# 'bim_strand',
# 'bim_site_coding_loc_base0',
# 'bim_site_loc_in_codon_base0',
# 'bim_original_codon',
# 'bim_nuc',
# 'bim_location_in_area',
# 'bim_site_area',
# 'bim_site_key',
# 'squ_bim_key',
# 'squ_bim_ortholog_sequences',
# 'squ_bim_same_target_aa',
# 'squ_bim_same_original_aa',
# 'squ_bim_same_position_in_codon',
# 'squ_bim_conserved',
# 'squ_bim_site_key']
        