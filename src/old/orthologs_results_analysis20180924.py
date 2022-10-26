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


"""
check if sequences in query_fasta are all contained within subject_fasta by id
and that they are all identical
"""
def compare_sequences(query_fasta, subject_fasta):
    
    cwd = os.getcwd()
    
    subbject_fasta_dict = {}
    for srecord in SeqIO.parse(open(subject_fasta, "r"), "fasta"):
        subject_fasta_dict.update({srecord.id:str(srecord.seq)})
        
    query_records_not_in_subject = []
    query_records_not_identical_to_subject = []
    with open(cwd + '/comparison.txt', 'w') as comparison_file:
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

"""
read a txt file with rows as:
animal|site_key
and create a dictionay of all keys per animal
"""
def keys_list_dictionary(file_path, file_name, split_by = '|'):
    
    animal_keys_dict = {}
    with open(file_path + file_name, "r") as f: 
        content = f.readlines()
        for line in content:
            key = line.split(split_by)
            animal = key[0]
            record_id = key[1].replace(' ','').replace('\n','')
            if key[0] in animal_keys_dict:
                animal_keys_dict[animal].append(record_id)
            else:
                animal_keys_dict.update({animal:[record_id]})
                        
    #print all list to excel sheets
    book = xlsxwriter.Workbook(file_path + 'animals_comps.xlsx')
    for key, val in animal_keys_dict.items():
        sheet = book.add_worksheet(key)
        for i, site in enumerate(val):
            sheet.write(i,0,site)
    book.close()
    
    return animal_keys_dict


#==============================================================================
"""
find regex (regex - pre compiled regex) in header
"""
def find_by_regex_in_header(header, regex):
    try:
        return regex.findall(header)[0]
    except:
        return 'unknown'  


def get_orthologs_data(input_path, match_file):
    
    data =[]
    
    #those columns repeat for two animals compared:
    #the raw data contain the animal and seq_id data as a continious string, seperated by the algorithm to two different strings
    #after thos columns(X2) there are two additional columns - alignment_identity and strand_orientation
    raw_columns_list = ['seq_id','protein','position','nuc_mm','aa_mm','editing_level','rna_coverage']
    
    #seperator for the animal|seq_id string
    animal_seqid_regex = re.compile(r'(.+?)\|([^\s]+)')
    n=0
    with open(input_path + match_file, "r") as f:
        content = f.readlines()    
        for i,line in enumerate(content):
            n+=1
            line_arr = line.split('\t')
            animal_seqid_1 = find_by_regex_in_header(line_arr[0],animal_seqid_regex)
            animal_seqid_2 = find_by_regex_in_header(line_arr[7],animal_seqid_regex)
            line_arr = [animal_seqid_1[1]] + line_arr[1:7] + [animal_seqid_2[1]] + line_arr[8:14] + [line_arr[14]]
            data.append(line_arr) 
        
        columns = [animal_seqid_1[0]+'_'+x for x in raw_columns_list] + [animal_seqid_2[0]+'_'+x for x in raw_columns_list] + [animal_seqid_1[0]+'_'+animal_seqid_2[0]+'_alignment_identity']
        
        orthologs_df = pd.DataFrame(data = data, columns = columns)
        orthologs_df = orthologs_df.apply(pd.to_numeric, errors='ignore')
        #create key columns for both animals (seq_id + site position)
        orthologs_df[animal_seqid_1[0]+'_key'] = orthologs_df.apply(lambda x: x[animal_seqid_1[0]+'_seq_id'] + '_' + str(x[animal_seqid_1[0]+'_position']), axis = 1)
        orthologs_df[animal_seqid_2[0]+'_key'] = orthologs_df.apply(lambda x: x[animal_seqid_2[0]+'_seq_id'] + '_' + str(x[animal_seqid_2[0]+'_position']), axis = 1)
        orthologs_df[animal_seqid_1[0]+'_'+animal_seqid_2[0]+'_recoding'] = orthologs_df.apply(lambda x: False if x[animal_seqid_1[0]+'_aa_mm'][0] == x[animal_seqid_1[0]+'_aa_mm'][1] and x[animal_seqid_2[0]+'_aa_mm'][0] == x[animal_seqid_2[0]+'_aa_mm'][1] else True, axis = 1)
        print(str(n) + ' matching sites in file')
     
    return orthologs_df, animal_seqid_1[0], animal_seqid_2[0]


"""
for QA purposes - retuen all conserved sites shared by results from find_all_matching_mm_fixed.pl
when executed with the parameters animal_1 animal_2
and when executed with the parameters animal_2 animal_1
"""
def check_intersection_of_pair_results_different_order(pair_df1,pair_df2,names = ['bob','lin']):
    
    l1 = list(pair_df1.apply(lambda x: x[names[0]+'_key'] + '_' + x[names[1]+'_key'], axis = 1)) 
    l2 = list(pair_df2.apply(lambda x: x[names[0]+'_key'] + '_' + x[names[1]+'_key'], axis = 1))
    return list(set(l1)&set(l2))

"""
shit elemnts in a list (arr) by n slots
"""
def shift(arr,n):
    return it.islice(it.cycle(arr),n,n+len(arr))


def scatter_orthologs_editing_levels(output_folder,el_1,el_2,animal_1,animal_2):

    pearson_corel_coef, p_val = stats.pearsonr(el_1, el_2)
    fig, ax = plt.subplots()
    plt.scatter(el_1,el_2, c = 'red', marker = 'D')
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
        
def plot_conserved_recoding_sites_fractions_per_pair(output_folder,all_dfs_dict,pairs_to_plot = []):

    vals = []
    labels = []
    
    if not len(pairs_to_plot):
        pairs_to_plot = list(all_dfs_dict.keys())
    
    for pair in pairs_to_plot:
        pair_df = all_dfs_dict[pair]
        labels.append(pair)
        vals.append(len(pair_df[pair_df[pair+'_recoding']])/float(len(all_dfs_dict[pair])))
     
    fig, ax = plt.subplots()
    y_pos = np.arange(len(pairs_to_plot))
    
    ax.barh(y_pos, vals, align='center', color='blue')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Fraction of conserved recoding sites')
    
    plt.show()
    plt.close()

def plot_conservations_sites_all_pairs(output_folder,all_dfs_dict, recoding = False, pairs_to_plot = []):
    
    if not len(pairs_to_plot):
        pairs_to_plot = list(all_dfs_dict.keys())
    
    mm_list = ['AC','AG','AT','CA','CG','CT','GA','GC','GT','TA','TC','TG']    
    colors = ['y','r','b','g','m','tan','silver','yellow','crimson','teal','k','orange','brown','limegreen']
    pairs_number = len(all_dfs_dict) 
    ind = np.arange(len(mm_list))   # the x locations for the groups
    width = 0.8/(len(all_dfs_dict)) # the width of the bars
    locations = np.arange(-pairs_number/2,pairs_number/2)
    axes_dict = {}
    animals_labels = []
    
    fig, ax = plt.subplots()
    
    for i, item in enumerate(all_dfs_dict.items()):  #iterate over all pairs of animals from dictionary
        if item[0] in pairs_to_plot:
            consereved_editing_types = []
            pair_df = item[1]
            animals = item[0].split('_')
            animals_labels.append(item[0].replace('_','-'))
#            print('pair ' + animals[0] + animals[1])
            for mm in mm_list: #count amount of each consereved mm type (if mm is the same for both animals in pair) 
                df_mm_view = pair_df[np.logical_and(pair_df[animals[0]+'_nuc_mm'] == mm,pair_df[animals[1]+'_nuc_mm'] == mm)].copy()
                if recoding:
                    if len(df_mm_view):    
                        df_mm_view['recoding'] = df_mm_view.apply(lambda x: False if x[animals[0]+'_aa_mm'][0] == x[animals[0]+'_aa_mm'][1] and x[animals[1]+'_aa_mm'][0] == x[animals[1]+'_aa_mm'][1] else True, axis = 1)
                        df_mm_view = df_mm_view[df_mm_view['recoding']]
                    else:
                        df_mm_view = []
                consereved_editing_types.append(len(df_mm_view)) 
            axes_dict[item[0]] = ax.bar(ind + width*locations[i], consereved_editing_types, width, color=colors[i])
#            print(consereved_editing_types)           
    
    ax.set_xticks(ind)
    ax.set_xticklabels(mm_list)
    ax.legend([x[0] for x in axes_dict.values()], animals_labels)
    
    if recoding:
        ax.set_title('Conserved Recoding Sites')
        fig.savefig(output_folder + 'number_of_conserves_recoding_sites_per_pair.png')
    else:
        ax.set_title('Conserved Sites')
        fig.savefig(output_folder + 'number_of_conserves_sites_per_pair.png')
    
    plt.close()
    
    

if __name__ == '__main__':
    
    input_path = 'C:/Users/user/Google_Drive/RNA_Editing/ortomcl_data/results/'
        
#    animals = ['sep','squ','oct','bim'] #animals strings as apear in files names (suppose to be identical to animals string within the file but not necessarily)
    animals = ['bob','lin','sep','squ','oct']
    animals_deq = deque(animals)
    raw_columns_list = ['seq_id','protein','position','nuc_mm','aa_mm','editing_level','rna_coverage']
    df_names_list = []
    all_dfs_dict = {}
    
    if not os.path.exists(input_path + 'analysis/'+'_'.join(animals)+'/'):
        os.makedirs(input_path + 'analysis/'+'_'.join(animals)+'/')
    output_folder = input_path + 'analysis/'+'_'.join(animals)+'/'    
    
    for i in range(len(animals)):
        animal_1 = animals[i]
        animals_deq.popleft()
#        print(animals_deq)
        for j in range(len(animals_deq)):
            animal_2 = animals_deq[j] 
            animals_names = [animal_1,animal_2]
            animals_sequenced = animal_1 + '_' + animal_2            
            print(animals_sequenced)
            match_file = 'matching_mms_'+animal_1+'_'+animal_2+'.txt'
            try:
                #create temporary df of matching sites from both animals (also return animals names as apear in file content)
                orthologs_df_temp, animal_1_name, animal_2_name = get_orthologs_data(input_path, match_file)
                #append pair df to pairs dictionary
                all_dfs_dict.update({animal_1 + '_' + animal_2:orthologs_df_temp})
                
                #if df_all_animals already contain animals, merge with temporary df. the result should yield a df with sites cosereved over all pairs of animals from from animals list
                if len(df_names_list):
                    #both animals already exist in df_all_animals - merge by both animals kyes and eventually take from temp df the alignment_identity column only
                    if all(a+'_key' in df_all_animals.columns for a in animals_names):
                        print('both animals are already in df')
                        columns_to_use = [animal_1_name+'_'+animal_2_name+'_alignment_identity'] + [animal_1_name + '_key',animal_2_name + '_key'] + [animal_1_name+'_'+animal_2_name+'_recoding']
                        df_all_animals = df_all_animals.merge(orthologs_df_temp[columns_to_use], on = [animal_1_name + '_key',animal_2_name + '_key'], how = 'inner')
                    #only animal1 exits in df_all_animals - merge by animal1 key and eventually take from temp df all animal2 columns (and alignment identity)
                    elif animal_1_name + '_key' in df_all_animals.columns:
                        print('animal_1 (' + animal_1_name + ') is already in df')
                        columns_to_use = [animal_2_name + '_' + x for x in raw_columns_list] + [animal_1_name + '_key',animal_2_name + '_key'] + [animal_1_name+'_'+animal_2_name+'_alignment_identity'] + [animal_1_name+'_'+animal_2_name+'_recoding']
                        df_all_animals = df_all_animals.merge(orthologs_df_temp[columns_to_use], on = animal_1_name + '_key', how = 'inner')
                    #same logic as last condition but now only animal2 is already exist in df_all_animals
                    elif animal_2_name + '_key' in df_all_animals.columns:
                        print('animal_2 (' + animal_2_name + ') is already in df')
                        columns_to_use = [animal_1_name + '_' + x for x in raw_columns_list] + [animal_1_name + '_key',animal_2_name + '_key'] + [animal_1_name+'_'+animal_2_name+'_alignment_identity'] + [animal_1_name+'_'+animal_2_name+'_recoding']
                        df_all_animals = df_all_animals.merge(orthologs_df_temp[columns_to_use], on = animal_2_name + '_key', how = 'inner')
                    if all(a+'_key' not in df_all_animals.columns for a in animals_names):
                        print('NO INTERSECTION OF ANIMALS')   
                else:
                    df_all_animals = orthologs_df_temp
                
                df_names_list.append([animal_1_name,animal_2_name])
                                
            except(FileNotFoundError):
                print('no file sequence ' + animal_1 + '_' + animal_2)
    
    
    #plot consereved sites for all animals (all, and recoding only)
    plot_conservations_sites_all_pairs(output_folder,all_dfs_dict,recoding = False)
    plot_conservations_sites_all_pairs(output_folder,all_dfs_dict,recoding = True)

    #for each pair of animals scatter editing level of sites shared by all animals    
    for pair in df_names_list:
        print(pair)
        df_all_animals_recoding = df_all_animals[df_all_animals[pair[0]+'_'+pair[1]+'_recoding']]
        print(len(df_all_animals_recoding))
        el_1 = df_all_animals_recoding[pair[0]+'_editing_level']
        el_2 = df_all_animals_recoding[pair[1]+'_editing_level']
        scatter_orthologs_editing_levels(output_folder,el_1,el_2,pair[0],pair[1])
        
    df_all_animals['all_recoding'] = df_all_animals.apply(lambda x: True if all([x[pair[0]+'_'+pair[1]+'_recoding'] for pair in df_names_list]) else False, axis = 1)
    df_all_animals_recoding = df_all_animals[df_all_animals['all_recoding']]
    
    #print all pairs dataframes (and sites conserved acrros all animals) to an excel file 
    excel_writer = pd.ExcelWriter(output_folder+'_'.join(animals)+'.xlsx',engine='xlsxwriter')
    for dataframe, sheet in zip(all_dfs_dict.values(),all_dfs_dict.keys()):
        dataframe.to_excel(excel_writer, sheet_name=sheet, startrow=0 , startcol=0, index=False) 
    df_all_animals.to_excel(excel_writer, sheet_name='coserved_across_all',startrow=0 , startcol=0, index=False)
    excel_writer.save()
